#' @title Calculating scale-free network based on correlation matrix
#' @description
#'  Calculate the adjacency matrix, topological matrix, and identify the optimal soft threshold for defining scale-free network.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param cordata a correlation matrix.
#' @param power.vec a vector of soft thresholding powers for which the scale free topology fit indices are to be calculated.
#' @param network.type network type. Allowed values are (unique abbreviations of) "unsigned" or "signed".
#' @param RsquaredCut desired minimum scale free topology fitting index R^2.
#' @export
lzq_calSFNetforCorMatrix <- function(
    cordata,
    power.vec = c(c(1:10), seq(from = 12, to = 20, by = 2)),
    network.type = "signed",
    RsquaredCut = 0.85) {
  tmpfun <- function(cordata, power, network.type = "signed") {
    if (network.type == "signed") {
      ADJ <- abs(0.5 + 0.5 * cordata)^power
    }
    if (network.type == "unsigned") {
      ADJ <- abs(cordata)^power
    }
    k <- apply(ADJ, 2, sum) - 1 # Connectivity
    cut1 <- cut(k, 10)
    binned.k <- tapply(k, cut1, mean) # Mean Connectivity for each Bin
    freq1 <- tapply(k, cut1, length) / length(k)
    xx <- as.vector(log10(binned.k))
    lm1 <- stats::lm(as.numeric(log10(freq1 + .000000001)) ~ xx)
    return(data.frame(
      Power = power,
      SFT.R.sq = summary(lm1)$r.squared,
      slope = lm1$coefficients[[2]],
      mean.k = mean(k)
    ))
  }
  pickSotfThres <- purrr::map_df(power.vec, \(x)tmpfun(cordata,
    power = x,
    network.type = network.type
  ))
  pick_power <- ifelse(max(pickSotfThres$SFT.R.sq) >= RsquaredCut,
    which(pickSotfThres$SFT.R.sq >= RsquaredCut)[1],
    pickSotfThres$Power[which.max(pickSotfThres$SFT.R.sq)]
  )

  if (network.type == "signed") {
    ADJ <- abs(0.5 + 0.5 * cordata)^pick_power
  }
  if (network.type == "unsigned") {
    ADJ <- abs(cordata)^pick_power
  }
  diag(ADJ) <- 0
  ADJ[is.na(ADJ)] <- 0

  numTOM <- stats::as.dist(ADJ %*% ADJ + ADJ)

  kk <- apply(ADJ, 2, sum) # Connectivity for each gene

  Dhelp1 <- matrix(kk, ncol = length(kk), nrow = length(kk))
  denomTOM <- pmin(stats::as.dist(Dhelp1), stats::as.dist(t(Dhelp1))) + stats::as.dist(1 - ADJ)
  TOMmatrix <- as.matrix(numTOM / denomTOM)

  return(list(
    pickSotfThres = pickSotfThres,
    pick_power = pick_power,
    Connectivity = kk,
    ADJmatrix = ADJ,
    TOMmatrix = TOMmatrix
  ))
}
