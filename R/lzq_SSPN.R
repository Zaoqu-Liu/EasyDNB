#' @title Sample-specific perturbation network
#' @description
#' Performing sample-specific perturbation network analysis.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param expr A expression dataframe with gene rows and sample columns.
#' @param ref.samples Samples for constructing the background reference network; Column ID or names.
#' @param cor.method specifies the method for correlation analysis.
#' @param p.adjust.method correction method, a character string. Can be abbreviated. c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param ppi Protein-protein interaction network; background network.
#' @param min.combined.score Minimum combined score for determining protein-protein interaction.
#' @param nCores The number of cores will be used.
#' @import future
#' @import magrittr
#' @export
lzq_SSPN <- function(
    expr,
    ref.samples,
    cor.method = "pearson",
    p.adjust.method = "BH",
    ppi = ppi_h,
    min.combined.score = 900,
    nCores = parallel::detectCores() - 10) {
  cat("Publish details: Liu X et al. Personalized characterization of diseases using sample-specific networks. Nucleic Acids Res. 2016 Dec 15;44(22):e164. doi: 10.1093/nar/gkw772\n")
  cat("\n")

  cat("+++ Performing sample-specific perturbation network analysis...\n")
  cat("\n")
  if (is.numeric(ref.samples)) {
    case.samples <- colnames(expr)[-ref.samples]
  } else {
    case.samples <- colnames(expr)[!colnames(expr) %in% ref.samples]
  }

  ppi2 <- ppi[ppi$G1 %in% rownames(expr) & ppi$G2 %in% rownames(expr) & ppi$combined_score >= min.combined.score, -3]
  expr <- expr[unique(c(ppi2$G1)), ]
  cat(paste0("+++ ", nrow(expr), " intersected genes between expression matrix and PPI network...\n"))
  cat("\n")

  ppi3 <- as.data.frame(t(apply(ppi2, 1, sort))) %>% dplyr::distinct(V1, V2, .keep_all = T)
  cat(paste0("+++ ", nrow(ppi3), " edges in PPI network...\n"))
  cat("\n")

  cat("+++ Calculating background network using reference samples...\n")
  cat("\n")
  refExpr <- expr[, ref.samples]
  refCor <- suppressWarnings(CorandPval(t(refExpr),
    method = cor.method,
    p.adjust.method = p.adjust.method
  ))$r
  ppi3$refCor <- purrr::map_vec(1:nrow(ppi3), ~ refCor[ppi3$V1[.x], ppi3$V2[.x]])

  tmp_fun <- function(id) {
    caseExpr <- cbind(refExpr, case = expr[, id])
    caseCor <- suppressWarnings(CorandPval(t(caseExpr), method = cor.method, p.adjust.method = p.adjust.method))$r
    caseCor2 <- purrr::map_vec(1:nrow(ppi3), ~ caseCor[ppi3$V1[.x], ppi3$V2[.x]])
    casePert <- caseCor2 - ppi3$refCor
    caseZ <- casePert / ((1 - ppi3$refCor^2) / (length(ref.samples) - 1))
    caseP <- purrr::map_vec(caseZ, function(x) {
      2 * stats::pnorm(-abs(x))
    })
    caseP[is.na(caseP)|is.nan(caseP)] <- 1
    caseFDR <- stats::p.adjust(caseP, p.adjust.method)
    return(data.frame(V1 = ppi3$V1, V2 = ppi3$V2, corPert = casePert, caseZ = caseZ, caseP = caseP, caseFDR = caseFDR))
  }

  cat("+++ Calculating sample-specific perturbation network for each case sample...\n")
  cat("\n")

  plan(multisession, workers = nCores)
  caseSSPN.list <- furrr::future_map(case.samples, ~ tmp_fun(.x), .progress = T)
  names(caseSSPN.list) <- case.samples

  case.SSPN.Pert.matrix <- as.data.frame(purrr::map_df(caseSSPN.list, ~ .x[, "corPert"]))
  case.SSPN.P.matrix <- as.data.frame(purrr::map_df(caseSSPN.list, ~ .x[, "caseP"]))
  case.SSPN.FDR.matrix <- as.data.frame(purrr::map_df(caseSSPN.list, ~ .x[, "caseFDR"]))
  rownames(case.SSPN.Pert.matrix) <- rownames(case.SSPN.P.matrix) <-
    rownames(case.SSPN.FDR.matrix) <- paste0(ppi3$V1, "_", ppi3$V2)

  cat("\n")
  cat("+++ Done!\n")

  return(list(
    case.SSPN.list = caseSSPN.list,
    refNet = ppi3,
    case.SSPN.Pert.matrix = case.SSPN.Pert.matrix,
    case.SSPN.P.matrix = case.SSPN.P.matrix,
    case.SSPN.FDR.matrix = case.SSPN.FDR.matrix,
    PPI.used = ppi3
  ))
}
