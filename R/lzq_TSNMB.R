#' @title Time-series network module biomarker analysis
#' @description
#' Performing time-series network module biomarker analysis.
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
#' @import magrittr
#' @export
lzq_TSNMB <- function(
    expr,
    state,
    state.levels,
    cor.method = "pearson",
    p.adjust.method = "BH",
    ppi = ppi_h,
    min.combined.score = 900,
    min.first.neighbor.size = 3,
    percent = T,
    top.n = 30,
    top.p = 0.05) {
  cat("+++ ThiS method was modified by Zaoqu Liu on the basis of sNMB!")
  cat("\n")

  cat("+++ Performing time-series network module biomarker analysis...\n")
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

  refSD <- apply(refExpr, 1, stats::sd)

  refLNSD <- refSD[names(ppil)]

  refCor <- suppressWarnings(CorandPval(t(refExpr),
    method = cor.method,
    p.adjust.method = p.adjust.method
  ))$r

  refNet <- purrr::map2(names(ppil), ppil, ~ data.frame(Gene = .y, PCC = refCor[.y, .x], SD = refSD[.y], row.names = NULL))
  names(refNet) <- names(ppil)

  state_case <- state[state$state != "ref", ]

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

  tmp_fun <- function(G) {
    caseExpr <- cbind(refExpr, case = expr[, state_case$ID[state_case$state == G]])
    caseSD <- apply(caseExpr, 1, stats::sd)
    caseLNSD <- caseSD[names(ppil)]
    DLNSD <- abs(caseLNSD - refLNSD)
    caseCor <- suppressWarnings(CorandPval(t(caseExpr), method = cor.method, p.adjust.method = p.adjust.method))$r
    caseNet <- purrr::map2(names(ppil), ppil, ~ data.frame(Gene = .y, PCC = caseCor[.y, .x], SD = caseSD[.y], row.names = NULL))
    names(caseNet) <- names(ppil)
    DNet <- purrr::map2(refNet, caseNet, \(x, y){
      data.frame(Gene = x$Gene, DPCC = abs(x$PCC - y$PCC), DSD = abs(x$SD - y$SD))
    })
    LI <- purrr::map2_vec(DNet, DLNSD, \(x, y){
      ((sum(x$DSD) + y) / (nrow(x) + 1)) * (sum(x$DPCC) / nrow(x))
    })
    LI <- data.frame(Gene = names(DLNSD), LI = LI, row.names = NULL) %>%
      dplyr::arrange(dplyr::desc(LI))
    GI <- sum(LI$LI[1:N]) / N
    return(list(DNet = DNet, LI = LI, GI = GI))
  }

  cat("+++ Calculating local sNMB score for each state...\n")
  cat("\n")
  state.res <- purrr::map(state.levels[-1], ~ tmp_fun(.x), .progress = T)
  names(state.res) <- state.levels[-1]
  state.LI.list <- purrr::map(state.res, ~ .x$LI)

  cat("\n")
  cat("+++ Calculating global sNMB score for each state...\n")
  cat("\n")
  state.GI <- data.frame(state = names(state.res), GI = purrr::map_vec(state.res, ~ .x$GI), row.names = NULL)

  state.DNet <- purrr::map(state.res, ~ .x$DNet)

  signaling.gene <- state.LI.list[[as.character(state.GI$state[which.max(state.GI$GI)])]] %>%
    dplyr::arrange(dplyr::desc(LI)) %>%
    dplyr::top_n(n = N, wt = LI) %>%
    .[, 1]

  cat(paste0("+++ ", state.GI$state[which.max(state.GI$GI)], " is the critical state...\n"))
  cat(paste0("+++ ", length(signaling.gene), " signaling genes involved in the critical state...\n"))
  cat("\n")
  cat("+++ Done!\n")

  return(list(
    state.GI = state.GI,
    state.LI.list = state.LI.list,
    signaling.gene = signaling.gene,
    state.DNet = state.DNet,
    refNet = refNet,
    refLNSD = refLNSD,
    Gene_module = ppil,
    PPI.used = ppi2
  ))
}
