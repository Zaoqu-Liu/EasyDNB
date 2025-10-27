#' @title Landscape DNB analysis
#' @description
#' Performing landscape dynamic network biomarker analysis.
#' @author Zaoqu Liu; E-mail: liuzaoqu@163.com
#' @param expr A expression dataframe with gene rows and sample columns.
#' @param state A time-series dataframe with two columns, the first is the sample names and the second is the group or time point information.
#' @param state.levels A vector for state sequence.
#' @param cor.method specifies the method for correlation analysis.
#' @param p.adjust.method correction method, a character string. Can be abbreviated. c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param variation.method specifies the method for calculating gene variation. sd or cv.
#' @param ppi Protein-protein interaction network; background network.
#' @param min.combined.score Minimum combined score for determining protein-protein interaction.
#' @param min.first.neighbor.size Minimum size of first order genes of a specific center gene.
#' @param min.second.neighbor.size Minimum size of second order genes of a specific center gene.
#' @param use.PCC.P.type Nominal P value (NP) or adjust P value (FDR) were selected to define the significance of sPCC.
#' @param use.PCC.P.cutoff A cutoff value for use.PCC.P.type.
#' @param percent Whether to use Percent to determine the number of DNB genes.
#' @param top.n Only Percent = F takes effect. Center genes with top (number) DNB score were defined as DNB genes.
#' @param top.p Only Percent = T takes effect. Center genes with top (percent) DNB score were defined as DNB genes.
#' @param nCores The number of cores will be used.
#' @import future
#' @import magrittr
#' @export
lzq_LDNB <- function(
    expr,
    state,
    state.levels,
    cor.method = "pearson",
    p.adjust.method = "BH",
    ppi = ppi_h,
    min.combined.score = 990,
    min.first.neighbor.size = 10,
    min.second.neighbor.size = 1,
    use.PCC.P.type = 'FDR',
    use.PCC.P.cutoff = 0.05,
    percent = F,
    top.n = 30,
    top.p = 0.05,
    nCores = parallel::detectCores() - 10
){
  cat("Publish details: Liu X et al. Detection for disease tipping points by landscape dynamic network biomarkers. Natl Sci Rev. 2019 Jul;6(4):775-785. doi: 10.1093/nsr/nwy162\n")
  cat("\n")
  cat("Publish details: Zhang C et al. Landscape dynamic network biomarker analysis reveals the tipping point of transcriptome reprogramming to prevent skin photodamage. J Mol Cell Biol. 2022 Jan 21;13(11):822-833. doi: 10.1093/jmcb/mjab060\n")
  cat("\n")

  cat("+++ Performing landscape dynamic network biomarker analysis...\n")
  cat("\n")
  if (ncol(state) != 2) {
    stop("Number of state columns must be two, the first is the sample names, and the second is the group or time point information!")
  }
  colnames(state) <- c("ID", "state")
  state$state <- factor(state$state, state.levels)
  state <- state[match(colnames(expr), state$ID), ]

  ref.samples <- state$ID[state$state=='ref']
  SSN <- lzq_SSPN(expr = expr,
                  ref.samples = ref.samples,
                  cor.method = cor.method,
                  p.adjust.method=p.adjust.method,
                  ppi=ppi,
                  min.combined.score=min.combined.score,
                  nCores=nCores)
  case.SSPN.list <- SSN$case.SSPN.list

  expr <- expr[unique(c(unique(case.SSPN.list[[1]]$V1),unique(case.SSPN.list[[1]]$V2))),]
  refExpr <- expr[,ref.samples]
  refMean <- rowMeans(refExpr)

  case.SSPN.list2 <- purrr::map(case.SSPN.list,\(x){
    if (use.PCC.P.type=='FDR') {
      x2 <- x[x$caseFDR<use.PCC.P.cutoff,]
    }else{
      x2 <- x[x$caseP<use.PCC.P.cutoff,]
    }
    x3 <- x2[,c(2,1,3:6)]
    colnames(x3)[1:2] <- c('V1','V2')
    return(rbind(x2,x3))
  })

  plan(multisession, workers = nCores)
  res.list <- furrr::future_map2(case.SSPN.list2,names(case.SSPN.list2),\(s,id){
    sppil <- purrr::map(unique(s$V1), \(V){
      return(s[s$V1 == V, 2])
    })
    names(sppil) <- unique(s$V1)
    retained <- purrr::map_vec(sppil, length) >= min.first.neighbor.size

    sppil <- sppil[retained] # first-order

    sExp <- expr[,id]; names(sExp) <- rownames(expr)
    DsExp <- abs(sExp-refMean)

    sED.in <- purrr::map2_vec(names(sppil),sppil,~mean(DsExp[c(.x,.y)]))
    sPCC.in <- purrr::map2_vec(names(sppil),sppil,~mean(abs(s[s$V1==.x&s$V2%in%.y,'corPert'])))
    sPCC.out <- purrr::map2_vec(names(sppil),sppil,\(x,y){
      s2 <- s[s$V1%in%y,]
      s2 <- s2[!s2$V2%in%c(x,y),]
      if(nrow(s2)>=min.second.neighbor.size){
        return(mean(abs(s2$corPert)))
      }else{
        return(NA)
      }
    })
    I <- sED.in * sPCC.in / sPCC.out
    retained <- !is.na(I)
    I <- I[retained]

    if (percent == F) {
      if (top.n > length(I)) {
        N <- length(I)
      }else{
        N <- top.n
      }
    } else {
      N <- ceiling(length(I) * top.p)
    }

    return(list(LI=data.frame(Gene=names(sppil)[retained],
                              sED.in=sED.in[retained],
                              sPCC.in=sPCC.in[retained],
                              sPCC.out=sPCC.out[retained],
                              LI=I),
                GI=mean(sort(I,decreasing = T)[1:N])))
  })

  case.LI.list <- purrr::map(res.list,~.x$LI)
  case.GI <- purrr::map2_df(res.list,names(res.list),~data.frame(ID=.y,GI=.x$GI))%>%
    dplyr::inner_join(x=state,y=.,by='ID')

  state.GI <- dplyr::group_by(case.GI,state)%>%
    dplyr::summarise(GI=mean(GI))

  # ?
  k <- case.LI.list[case.GI$ID[case.GI$state==as.character(state.GI$state[which.max(state.GI$GI)])]]
  candidates <- Reduce(intersect,purrr::map(k,~.x$Gene))
  k2 <- purrr::map(k,~.x[match(candidates,.x$Gene),c('Gene','LI')])
  k3 <- rowMeans(purrr::map_df(k2,~.x$LI))
  k4 <- data.frame(Gene=candidates,LI=k3)%>%
    dplyr::arrange(dplyr::desc(LI))

  if (percent == F) {
    if (top.n > nrow(k4)) {
      N <- nrow(k4)
    }else{
      N <- top.n
    }
  } else {
    N <- ceiling(nrow(expr) * top.p)
    if (N>nrow(k4)) {N <- nrow(k4)}
  }
  DNB.genes <- k4$Gene[1:N]

  cat("\n")
  cat(paste0("+++ ", state.GI$state[which.max(state.GI$GI)], " is the critical state...\n"))
  cat(paste0("+++ ", length(DNB.genes), " genes (DNB module) involved in the critical state...\n"))
  cat("\n")
  cat("+++ Done!\n")
  future::plan(future::sequential)
  return(list(
    state.GI = state.GI,
    DNB.genes = DNB.genes,
    Gene.LI = k4,
    case.GI = case.GI,
    case.LI.list = case.LI.list,
    SSPN = SSN
  ))
}
