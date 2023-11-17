#' @title Conventional DNB analysis
#' @description
#' Performing conventional dynamic network biomarker analysis.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param expr A expression dataframe with gene rows and sample columns.
#' @param state A time-series dataframe with two columns, the first is the sample names and the second is the group or time point information.
#' @param state.levels A vector for state sequence.
#' @param cor.method specifies the method for correlation analysis.
#' @param p.adjust.method correction method, a character string. Can be abbreviated. c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param variation.method specifies the method for calculating gene variation. sd or cv.
#' @param min.size Minimum gene number of gene modules.
#' @param max.size Maximum gene number of gene modules.
#' @param AddModuleSize Whether to consider gene module size when calculating DNB score.
#' @export
lzq_cDNB <- function(
    expr,
    state,
    state.levels,
    cor.method = "pearson",
    p.adjust.method = "BH",
    variation.method = "sd",
    min.size = 10,
    max.size = 2000,
    AddModuleSize = F) {
  cat("Publish details: Chen L et al. Detecting early-warning signals for sudden deterioration of complex diseases by dynamical network biomarkers. Sci Rep. 2012;2:342. doi: 10.1038/srep00342\n")
  cat("\n")

  cat("+++ Performing conventional dynamic network biomarker analysis...\n")
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
  corl <- suppressWarnings(purrr::map(ddl, ~ CorandPval(t(.x),
    method = cor.method,
    p.adjust.method = p.adjust.method
  )))

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


  cat("+++ Performing Hierarchical clustering for all genes ...\n")
  cat("\n")

  max.size <- ifelse(max.size >= length(vl[[1]]), length(vl[[1]]) * 0.8, max.size)
  modulel <- purrr::map(corl, \(x){
    dend <- stats::as.dendrogram(stats::hclust(stats::as.dist(1 - x$r)))
    modules <- dendextend::partition_leaves(dend)
    member_nums <- sapply(modules, length)
    modules <- modules[member_nums >= min.size & member_nums <= max.size]
    return(modules)
  })

  cat("+++ Number of gene modules detected in different groups or times:\n")
  cat(paste0(paste0("+++ ", names(modulel), " = ", sapply(modulel, length), collapse = "\n")), "\n")
  cat("\n")
  cat("+++ Calculating module scores for each group or time point...\n")
  cat("\n")

  score_l <- purrr::map2(data, modulel, \(x, y){
    correlation <- x$r
    v <- x$v
    mol <- y
    scores <- purrr::map_df(mol, \(genes_in){
      genes_out <- setdiff(names(v), genes_in)
      cor_in <- mean(stats::as.dist(correlation[genes_in, genes_in]))
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
      Module = 1:length(mol), Size = sapply(mol, length), CI = scores$s,
      V_in = scores$v_in, R_in = scores$r_in, R_out = scores$r_out
    ))
  })

  Candidate <- data.frame(
    State = names(score_l),
    Module = purrr::map_vec(score_l, ~ which.max(.x$CI)),
    CI = purrr::map_vec(score_l, ~ max(.x$CI)),
    V_in = purrr::map_vec(score_l, ~ .x$V_in[which.max(.x$CI)]),
    R_in = purrr::map_vec(score_l, ~ .x$R_in[which.max(.x$CI)]),
    R_out = purrr::map_vec(score_l, ~ .x$R_out[which.max(.x$CI)]),
    row.names = NULL
  )

  DNB.genes <- modulel[[Candidate$State[which.max(Candidate$CI)]]][[Candidate$Module[which.max(Candidate$CI)]]]

  DNB.score <- purrr::map_df(data, \(x){
    correlation <- x$r
    v <- x$v
    genes_in <- DNB.genes
    genes_out <- setdiff(names(v), genes_in)
    cor_in <- mean(stats::as.dist(correlation[genes_in, genes_in]))
    cor_out <- mean(correlation[genes_in, genes_out])
    v_in <- mean(v[genes_in])
    if (AddModuleSize) {
      score <- sqrt(length(genes_in)) * v_in * cor_in / cor_out
    } else {
      score <- v_in * cor_in / cor_out
    }
    return(data.frame(
      CI = score, V_in = v_in,
      R_in = cor_in, R_out = cor_out
    ))
  })
  DNB.score <- cbind(State = state.levels, DNB.score)

  cat(paste0("+++ ", Candidate$State[which.max(Candidate$CI)], " is the critical state...\n"))
  cat(paste0("+++ ", length(DNB.genes), " genes (DNB module) involved in the critical state...\n"))
  cat("\n")
  cat("+++ Done!\n")

  return(list(
    DNB.score = DNB.score, DNB.genes = DNB.genes,
    CI_all = score_l, Gene_module = modulel,
    Candidate = Candidate, Cor = corl, V = vl
  ))
}
