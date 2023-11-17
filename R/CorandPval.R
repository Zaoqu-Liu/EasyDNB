#' @title Calculation of correlations and associated p-values
#' @description
#' A faster, one-step calculation of Student correlation p-values for multiple correlations, properly taking into account the actual number of observations.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param x a vector or a matrix
#' @param y a vector or a matrix. If NULL, the correlation of columns of x will be calculated.
#' @param use determines handling of missing data. See cor() for details.
#' @param alternative specifies the alternative hypothesis and must be (a unique abbreviation of) one of "two.sided", "greater" or "less". the initial letter. "greater" corresponds to positive association, "less" to negative association.
#' @param p.adjust.method correction method, a character string. Can be abbreviated. c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param ... other arguments to the function cor().
#' @details
#' The function calculates correlations of a matrix or of two matrices and the corresponding Student p-values. The output is not as full-featured as cor.test(), but can work with matrices as input.
#' @export
CorandPval <- function(x, y = NULL, use = "everything",
                       r.na.value = 0, p.na.value = 1,
                       p.adjust.method = "BH", alternative = c("two.sided", "less", "greater"),
                       ...) {
  ia <- match.arg(alternative)
  cor <- cor(x, y, use = use, ...)
  x <- as.matrix(x)
  finMat <- !is.na(x)
  if (is.null(y)) {
    np <- t(finMat) %*% finMat
  } else {
    y <- as.matrix(y)
    np <- t(finMat) %*% (!is.na(y))
  }
  Z <- 0.5 * log((1 + cor) / (1 - cor)) * sqrt(np - 2)
  if (ia == "two.sided") {
    T <- sqrt(np - 2) * abs(cor) / sqrt(1 - cor^2)
    p <- 2 * stats::pt(T, np - 2, lower.tail = FALSE)
  } else if (ia == "less") {
    T <- sqrt(np - 2) * cor / sqrt(1 - cor^2)
    p <- stats::pt(T, np - 2, lower.tail = TRUE)
  } else if (ia == "greater") {
    T <- sqrt(np - 2) * cor / sqrt(1 - cor^2)
    p <- stats::pt(T, np - 2, lower.tail = FALSE)
  }

  cor[is.na(cor)] <- r.na.value
  p[is.na(p)] <- p.na.value

  fdr <- apply(p, 2, \(x)stats::p.adjust(x, p.adjust.method))

  list(r = cor, p = p, fdr = fdr)
}
