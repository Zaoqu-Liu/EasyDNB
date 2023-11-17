#' @title Single-sample landscape entropy analysis
#' @description
#' Performing single-sample landscape entropy analysis.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param expr A expression dataframe with gene rows and sample columns.
#' @param state A time-series dataframe with two columns, the first is the sample names and the second is the group or time point information.
#' @param state.levels A vector for state sequence.
#' @param cor.method specifies the method for correlation analysis.
#' @param p.adjust.method correction method, a character string. Can be abbreviated. c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param min.first.neighbor.size Minimum size of first order genes of a specific center gene.
#' @param ppi Protein-protein interaction network; background network.
#' @param min.combined.score Minimum combined score for determining protein-protein interaction.
#' @param percent Whether to use Percent to determine the number of DNB genes.
#' @param top.n Only Percent = F takes effect. Center genes with top (number) DNB score were defined as DNB genes.
#' @param top.p Only Percent = T takes effect. Center genes with top (percent) DNB score were defined as DNB genes.
#' @param nCores The number of cores will be used.
#' @import future
#' @import magrittr
#' @export
lzq_SLE <- function(
    expr,
    state,
    state.levels,
    cor.method = "pearson",
    p.adjust.method = "BH",
    ppi = ppi_h,
    min.combined.score = 900,
    min.first.neighbor.size = 1,
    percent = T,
    top.n = 30,
    top.p = 0.05,
    nCores = parallel::detectCores() - 10) {
  # reference: Single-sample landscape entropy reveals the imminent phase transition during disease progression

  cat("Publish details: Liu R et al. Single-sample landscape entropy reveals the imminent phase transition during disease progression. Bioinformatics. 2020 Mar 1;36(5):1522-1532. doi: 10.1093/bioinformatics/btz758\n")
  cat("\n")

  cat("+++ Performing single-sample landscape entropy analysis...\n")
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

  cat(paste0("+++ ", length(ppil), " centre genes detected in this dataset...\n"))
  cat("\n")
  cat("+++ Calculating background network using reference samples...\n")
  cat("\n")

  refExpr <- expr[, state$ID[state$state == "ref"]]
  refCor <- suppressWarnings(CorandPval(t(refExpr),
    method = cor.method,
    p.adjust.method = p.adjust.method
  ))$r
  refLE <- purrr::map2_df(ppil, names(ppil), \(x, y){
    E2 <- abs(refCor[x, y]) / sum(abs(refCor[x, y])) # Equation 2
    E2 <- E2[E2 != 0]
    E1 <- -sum(E2 * log(E2)) / log(length(E2)) # Equation 1
    RSD <- stats::sd(as.numeric(refExpr[y, ]))
    return(data.frame(Gene = y, RLE = E1, RSD = RSD))
  }) %>%
    dplyr::filter(!is.na(RLE), RSD > 0)

  ppil <- ppil[refLE$Gene]

  state_case <- state[state$state != "ref", ]

  tmp_fun <- function(id) {
    caseExpr <- cbind(refExpr, case = expr[, id])
    caseCor <- suppressWarnings(CorandPval(t(caseExpr), method = cor.method, p.adjust.method = p.adjust.method))$r
    caseLE <- purrr::map2_df(ppil, names(ppil), \(x, y){
      E2 <- abs(caseCor[x, y]) / sum(abs(caseCor[x, y])) # Equation 2
      E2 <- E2[E2 != 0]
      E3 <- -sum(E2 * log(E2)) / log(length(E2)) # Equation 3
      CSD <- stats::sd(as.numeric(caseExpr[y, ]))
      DSD <- abs(CSD - refLE$RSD[refLE$Gene == y]) # Equation 5
      DE <- abs(E3 - refLE$RLE[refLE$Gene == y])
      DH <- DE * DSD # Equation 4; The differential standard deviation is regarded as the weight coefficient.
      return(data.frame(Gene = y, CLE = E3, CSD = CSD, DLE = DE, DSD = DSD, LSLE = DH))
    })
  }

  cat("+++ Calculating local SLE score for each case sample...\n")
  cat("\n")
  plan(multisession, workers = nCores)
  caseLE.list <- furrr::future_map(state_case$ID, ~ tmp_fun(.x), .progress = T)
  names(caseLE.list) <- state_case$ID

  case.GLSE <- purrr::map2_df(caseLE.list, names(caseLE.list), \(x, y){
    data.frame(ID = y, GLSE = mean(x$LSLE))
  })

  cat("\n")
  cat("+++ Calculating global SLE score for each case sample...\n")
  cat("\n")
  case.GLSE <- merge(state_case, case.GLSE, by = 1)
  caseLE.integrated <- purrr::map(state.levels[-1], \(tt){
    l <- caseLE.list[case.GLSE$ID[case.GLSE$state == tt]]
    l2 <- rowMeans(purrr::map_df(l, ~ .x$LSLE)) # Equation 6
    data.frame(Gene = names(ppil), LSLE = l2)
  })
  names(caseLE.integrated) <- state.levels[-1]

  cat("+++ Calculating global SLE score for each state...\n")
  cat("\n")
  state.GLSE <- dplyr::group_by(case.GLSE, state) %>%
    dplyr::summarise(GLSE = mean(GLSE))

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

  SLE.gene <- caseLE.integrated[[as.character(state.GLSE$state[which.max(state.GLSE$GLSE)])]] %>%
    dplyr::arrange(dplyr::desc(LSLE)) %>%
    dplyr::top_n(n = N, wt = LSLE) %>%
    .[, 1]

  cat(paste0("+++ ", state.GLSE$state[which.max(state.GLSE$GLSE)], " is the critical state...\n"))
  cat(paste0("+++ ", length(SLE.gene), " genes (DNB module) involved in the critical state...\n"))
  cat("\n")
  cat("+++ Done!\n")

  return(list(
    GSLE.score = state.GLSE,
    SLE.gene = SLE.gene,
    stateLE.list = caseLE.integrated,
    caseLE.list = caseLE.list,
    refLE = refLE,
    Gene_module = ppil,
    PPI.used = ppi2
  ))
}
