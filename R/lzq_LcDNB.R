#' @title Landscape conventional DNB analysis
#' @description
#' Performing landscape conventional dynamic network biomarker analysis.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param expr A expression dataframe with gene rows and sample columns.
#' @param state A time-series dataframe with two columns, the first is the sample names and the second is the group or time point information.
#' @param state.levels A vector for state sequence.
#' @param cor.method specifies the method for correlation analysis.
#' @param p.adjust.method correction method, a character string. Can be abbreviated. c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param variation.method specifies the method for calculating gene variation. sd or cv.
#' @param min.first.neighbor.size Minimum size of first order genes of a specific center gene.
#' @param min.second.neighbor.size Minimum size of second order genes of a specific center gene.
#' @param ppi Protein-protein interaction network; background network.
#' @param min.combined.score Minimum combined score for determining protein-protein interaction.
#' @param percent Whether to use Percent to determine the number of DNB genes.
#' @param top.n Only Percent = F takes effect. Center genes with top (number) DNB score were defined as DNB genes.
#' @param top.p Only Percent = T takes effect. Center genes with top (percent) DNB score were defined as DNB genes.
#' @param AddModuleSize Whether to consider gene module size when calculating DNB score.
#' @export
lzq_LcDNB <- function(
    expr,
    state,
    state.levels,
    cor.method = "pearson",
    p.adjust.method = "BH",
    variation.method = "sd",
    min.first.neighbor.size = 3,
    min.second.neighbor.size = 1,
    ppi = ppi_h,
    min.combined.score = 900,
    percent = T,
    top.n = 30,
    top.p = 0.05,
    AddModuleSize = F) {
  cat("+++ ThiS method was modified by Zaoqu Liu on the basis of L-DNB!")
  cat("\n")

  cat("+++ Performing landscape conventional dynamic network biomarker analysis...\n")
  cat("\n")
  if (ncol(state) != 2) {
    stop("Number of state columns must be two, the first is the sample names, and the second is the group or time point information!")
  }
  colnames(state) <- c("ID", "state")
  state$state <- factor(state$state, state.levels)
  state <- state[match(colnames(expr), state$ID), ]
  ddl <- purrr::map(state.levels, \(x){
    expr[, state$ID[state$state == x]]
  })
  names(ddl) <- state.levels

  cat("+++ Calculating gene correlations for each group or time point...\n")
  cat("\n")
  corl <- suppressWarnings(purrr::map(ddl, ~ CorandPval(t(.x), method = cor.method, p.adjust.method = p.adjust.method)))

  cat("+++ Calculating gene varaitions for each group or time point...\n")
  cat("\n")
  if (variation.method == "sd") {
    vl <- purrr::map(ddl, ~ apply(.x, 1, sd))
  }
  if (variation.method == "cv") {
    vl <- purrr::map(ddl, ~ apply(.x, 1, \(a){
      sd(a) / mean(a)
    }))
  }

  data <- purrr::map2(corl, vl, \(x, y){
    list(r = x$r, v = y)
  })

  ppi2 <- ppi[ppi$G1 %in% rownames(expr) & ppi$G2 %in% rownames(expr) & ppi$combined_score >= min.combined.score, ]

  expr <- expr[unique(c(ppi2$G1)), ]
  cat(paste0("+++ ", nrow(expr), " intersected genes between expression matrix and PPI network...\n"))
  cat("\n")

  ppi3 <- ppi2[ppi2$G1 %in% names(table(ppi2$G1))[table(ppi2$G1) >= min.first.neighbor.size], ]
  ppil <- purrr::map(unique(ppi3$G1), \(x){
    g <- ppi2[ppi2$G1 == x, ]
    return(g$G2)
  })
  names(ppil) <- unique(ppi3$G1)

  second_detials <- purrr::map2(ppil, names(ppil), \(x, y){
    second_genes_df <- ppi2[unique(ppi2$G1) %in% c(x, y), ]
    second_genes <- unique(c(second_genes_df$G1, second_genes_df$G2))
    second_genes <- second_genes[!second_genes %in% c(x, y)]
  })
  retained <- purrr::map_vec(second_detials, length) >= min.second.neighbor.size

  ppil <- ppil[retained] # first-order
  second_detials <- second_detials[retained] # second-order

  cat(paste0("+++ ", length(ppil), " centre genes detected in this dataset...\n"))
  cat("\n")
  cat("+++ Calculating module scores for each group or time point...\n")
  cat("\n")

  score_l <- purrr::map(data, \(x){
    correlation <- x$r
    v <- x$v
    mol <- ppil
    scores <- purrr::map2_df(mol, names(mol), \(genes_in, core_gene){
      genes_out <- second_detials[[core_gene]]

      cor_in <- mean(correlation[core_gene, genes_in])
      cor_out <- mean(correlation[genes_in, genes_out])
      v_in <- mean(v[genes_in])

      if (AddModuleSize) {
        score <- sqrt(length(genes_in)) * v_in * cor_in / cor_out
      } else {
        score <- v_in * cor_in / cor_out
      }
      return(data.frame(s = score, v_in = v_in, r_in = cor_in, r_out = cor_out))
    })
    return(data.frame(
      Module = names(mol), Rank = -rank(scores$s) + length(mol) + 1,
      Size = sapply(mol, length), CI = scores$s,
      V_in = scores$v_in, R_in = scores$r_in, R_out = scores$r_out,
      row.names = NULL
    ))
  })

  if (percent == F) {
    if (top.n > length(ppil)) {
      stop(paste0(
        "\nThe number of candidata core genes is ",
        length(ppil), ", which is larger than the top.n you inputted!\n",
        "Please change another number!"
      ))
    }
    N <- top.n
  } else {
    N <- ceiling(length(ppil) * top.p)
  }

  Candidate <- data.frame(
    State = names(score_l),
    CI = purrr::map_vec(score_l, \(x){
      mean(x$CI[x$Rank %in% 1:N])
    }),
    row.names = NULL
  )

  Candidate_DNB.genes_list <- purrr::map(score_l, \(x){
    x$Module[x$Rank %in% 1:N]
  })
  DNB.genes <- Candidate_DNB.genes_list[[which.max(Candidate$CI)]]

  DNB.score <- purrr::map_vec(score_l, \(x){
    x2 <- x[x$Module %in% DNB.genes, ]
    mean(x2$CI)
  })
  DNB.score <- data.frame(State = state.levels, DNB.score = DNB.score)

  cat(paste0("+++ ", Candidate$State[which.max(Candidate$CI)], " is the critical state...\n"))
  cat(paste0("+++ ", length(DNB.genes), " genes (DNB module) involved in the critical state...\n"))
  cat("\n")
  cat("+++ Done!\n")

  return(list(
    DNB.score = DNB.score, DNB.genes = DNB.genes,
    CI_all = score_l, Gene_module = ppil,
    Candidate = Candidate, Cor = corl, V = vl, PPI.used = ppi2,
    first.order.genes = ppil, second.order.genes = second_detials
  ))
}
