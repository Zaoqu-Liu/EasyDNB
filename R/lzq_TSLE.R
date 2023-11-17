#' @title Time-series landscape entropy analysis
#' @description
#' Performing time-series landscape entropy analysis.
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
lzq_TSLE <- function(
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
  cat("+++ ThiS method was modified by Zaoqu Liu on the basis of SLE!")
  cat("\n")

  cat("+++ Performing time-series landscape entropy analysis...\n")
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
  cat("+++ Calculating background network for reference samples...\n")
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

  tmp_fun <- function(G) {
    caseExpr <- cbind(refExpr, case = expr[, state_case$ID[state_case$state == G]])
    caseCor <- suppressWarnings(CorandPval(t(caseExpr), method = cor.method, p.adjust.method = p.adjust.method))$r
    caseLE <- purrr::map2_df(ppil, names(ppil), \(x, y){
      E2 <- abs(caseCor[x, y]) / sum(abs(caseCor[x, y])) # Equation 2
      E2 <- E2[E2 != 0]
      E3 <- -sum(E2 * log(E2)) / log(length(E2)) # Equation 3
      CSD <- stats::sd(as.numeric(caseExpr[y, ]))
      DSD <- abs(CSD - refLE$RSD[refLE$Gene == y]) # Equation 5
      DE <- abs(E3 - refLE$RLE[refLE$Gene == y])
      DH <- DE * DSD # Equation 4; The differential standard deviation DSDðk; tÞ is regarded as the weight coefficient.
      return(data.frame(Gene = y, CLE = E3, CSD = CSD, DLE = DE, DSD = DSD, LSLE = DH))
    })
  }

  cat("+++ Calculating local SLE score for each state...\n")
  cat("\n")
  caseLE.list <- purrr::map(state.levels[-1], ~ tmp_fun(.x), .progress = T)
  names(caseLE.list) <- state.levels[-1]

  cat("+++ Calculating global SLE score for each state...\n")
  cat("\n")
  case.GLSE <- purrr::map2_df(caseLE.list, names(caseLE.list), \(x, y){
    data.frame(state = y, GLSE = mean(x$LSLE))
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

  SLE.gene <- caseLE.list[[as.character(case.GLSE$state[which.max(case.GLSE$GLSE)])]] %>%
    dplyr::arrange(dplyr::desc(LSLE)) %>%
    dplyr::top_n(n = N, wt = LSLE) %>%
    .[, 1]

  cat(paste0("+++ ", case.GLSE$state[which.max(case.GLSE$GLSE)], " is the critical state...\n"))
  cat(paste0("+++ ", length(SLE.gene), " genes (DNB module) involved in the critical state...\n"))
  cat("\n")
  cat("+++ Done!\n")

  return(list(
    GSLE.score = case.GLSE,
    SLE.gene = SLE.gene,
    stateLE.list = caseLE.list,
    refLE = refLE,
    Gene_module = ppil,
    PPI.used = ppi2
  ))
}
