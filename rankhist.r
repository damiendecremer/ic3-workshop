#########################################
#                                       # 
# RANK HISTOGRAM FOR ENSEMBLE FORECASTS #
#                                       #
#########################################
rankhist <- function(ens, ver) {
#
# calculate the rank histogram of a collection of ensemble forecasts
#
# Usage: 
#   rh <- rankhist(ens=ens, ver=ver)
#
# Arguments:
#
#   ens ... N*K matrix, rows are the ensemble forecasts
#   ver ... N vector of corresponding verifications
#
# Return value:
#   a vector of verification rank frequencies
#
# Author:
#   Stefan Siegert 
#   s.siegert@exeter.ac.uk 
#   October 2013
#
# Example:
#   ens <- matrix(rnorm(500), 100, 5)
#   ver <- rnorm(100)
#   rh <- rankhist(ens, ver)
#
# References: 
#   Talagrand (1997)
#   Hammill (2001)
#
#
  N <- dim(ens)[1]
  K <- dim(ens)[2]
  stopifnot(N == length(ver))
  ranks <- apply(cbind(ver, ens), 1, rank, ties.method="random")[1, ]
  rank.hist <- hist(ranks, breaks=seq(0.5, K+1.5, 1), plot=FALSE)$counts
  return(rank.hist)
}



#####################################
#                                   #
# RANK HISTOGRAM SIGNIFICANCE TESTS #
#                                   #
#####################################
rankhist.tests <- function(rank.hist) {
#
# Conduct a series of significance tests for flatness of the rank histogram
#
# Usage: rankhist.tests(rank.hist=rh)
#
# Arguments:
#   * rank.hist ... a vector of rank counts (see function `rankhist()`
#
# Return value:
#   * a list with the following elements
#     + pearson ... a dataframe with the pearson chi^2 test statistic 
#                   and its p-value under the chi^2 distribution with 
#                   K degrees of freedom
#     + jolliffe.primo ... a data frame, each column has the test 
#                          statistic and p-value of tests for a specific 
#                          alternative, namely for a sloped and convex
#                          rank histogram
#
# Author: 
#
#    Stefan Siegert 
#    s.siegert@exeter.ac.uk 
#    October 2013
#
# Example:
# 
#   ens <- matrix(rnorm(500), 100, 5)
#   ver <- rnorm(100)
#   rh <- rankhist(ens, ver)
#   rankhist.tests(rank.hist = rh)
#
# References: Pearson 1900
#             Jolliffe & Primo 2008
#             Broecker 2008
#
  o.i <- rank.hist
  N <- sum(o.i)
  J <- length(o.i) 
  i <- 1:J
  e.i <- N / J
  x.i <- (o.i - e.i) / sqrt(e.i)
  # pearson chi^2 test 
  X2 <- sum(x.i * x.i)
  p.chisq <- pchisq(X2, df=J-1, lower.tail=FALSE)
  # jolliffe-primo: 
  # linear contrast
  a <- 2 * sqrt(3 / (J^3 - J))
  b <- -(sqrt(3) * J + sqrt(3)) / sqrt(J * (J + 1) * (J - 1))
  x.lin <- a*i+b
  X2.lin <- sum(x.i * x.lin)^2 # should have chi^2(df=1)
  # squared contrast
  a <- 6 * sqrt(5 / (J^5 - 5 * J^3 + 4 * J))
  b <- -1 / 2 * (sqrt(5) * J^2 - sqrt(5)) / 
       (sqrt((J - 2) * (J - 1) * J * (J + 1) * (J + 2)))
  x.u <- a * (i - (J + 1) / 2)^2 + b
  X2.u <- sum(x.i * x.u)^2 # should have chisq(df=1)
  # return 
  pearson <- data.frame(test.statistic=X2, p.value=p.chisq)
  ret.df <- data.frame(pearson.chi2=c(X2, p.chisq), jp.slope=c(X2.lin, pchisq(X2.lin, df=1, lower.tail=FALSE)), jp.convex=c(X2.u, pchisq(X2.u, df=1, lower.tail=FALSE)))
  rownames(ret.df) <- c("test.statistic", "p.value")
  ret.df
}



#######################################
#                                     #
# RANK HIST COMPARISON BETWEEN TWO    #
#  ENSEMBLE FORECASTING SYSTEMS       #
#                                     #
#######################################
AnalyzeRankhistDifference <- function(ens, ens.ref, ver, n.boot=500) {
#
# Compare two ensembles in terms of their rank histograms. Statistics 
# that quantify the deviation from flatness of each rank histogram are 
# calculated. Their differences are bootstrapped to estimate their 
# sampling distributions. 
#
# Note: this function is still experimental and not recommended for public
# release
#
# Usage: 
#   AnalyzeRankhistDifference(ens=ens, enr.ref=ens.ref, ver=ver)
#
# Arguments:
#   * ens ... a N*K matrix; the ensemble that is being analyzed
#   * ens.ref ... a N*K.ref matrix; a reference ensemble to which `ens` is compared
#   * ver ... vector of length N; verifications to which the two ensemble forecasts refer
#   * n.boot ... number of repetitions of the resampling protocols
#
# Return value:
#   * a data.frame with three columns corresponding to the pearson chi^2
#     statistic, and two jolliffe-primo statistics that are sensitive to sloped
#     and convex rank histograms, respectively; the rows of the data.frame are as
#     follows:
#      + score.diff ... the observed difference between the two ensembles
#      + p.value ... the upper tail probability of observing `score.diff` under
#                    the bootstrapped distribution of no difference between
#                    the two ensembles; see `Details` for how this distribution
#                    is estimated
#      + Q0.01, etc ... a number of quantiles of the sampling distribution of 
#                       `score.diff`; see `Details` for how this distribution
#                        is estimated
#
# Details:
#   * the null-distribution of no difference is estimated by repeatedly
#     shuffling the ensemble members between the two models, and calculating the
#     three score differences between the two new ensembles; by shuffling the
#     ensemble members around, two artificial ensembles are constructed none of
#     which is superior to the other
#   * the sampling distribution of the score difference is estimated by
#     resampling N times with replacement from the time indices and calculating
#     the score differences between the ensembles at these time indices; possible
#     differences between the two ensembles are preserved; possibly not all
#     sources of randomness are accounted for which would make the resulting
#     confidence intervals too narrow
#
#
# Author: 
#
#    Stefan Siegert 
#    s.siegert@exeter.ac.uk 
#    October 2013
#
# Example:
# 
#   ens <- matrix(rnorm(500), 100, 5)
#   ver <- rnorm(100)
#   rh <- rankhist(ens, ver)
#   rankhist.tests(rank.hist = rh)
#
# References: Pearson 1900
#             Jolliffe & Primo 2008
#             Broecker 2008
#
# ens     ... ensemble
# ens.ref ... reference ensemble forecast
# ver     ... verification
# size    ... size of the one-sided statistical test
#
################################
  N <- nrow(ens)
  K <- ncol(ens)
  K.ref <- ncol(ens.ref)
  J <- K + 1
  J.ref <- K.ref + 1
  i <- 1:J
  i.ref <- 1:J.ref
  e.i <- N / J
  e.i.ref <- N / J.ref
  # linear contrasts
  a.lin <- 2 * sqrt(3 / (J^3 - J))
  b.lin <- -(sqrt(3) * (J + 1)) / sqrt(J * (J + 1) * (J - 1))
  c.lin <- a.lin * i + b.lin
  a.lin.ref <- 2 * sqrt(3 / (J.ref^3 - J.ref))
  b.lin.ref <- -(sqrt(3) * (J.ref + 1)) / sqrt(J.ref * (J.ref + 1) * (J.ref - 1))
  c.lin.ref <- a.lin.ref * i.ref + b.lin.ref
  # U-shaped contrasts
  a.u <- 6 * sqrt(5 / (J^5 - 5 * J^3 + 4 * J))
  b.u <- -1 / 2 * (sqrt(5) * J^2 - sqrt(5)) / 
       (sqrt((J - 2) * (J - 1) * J * (J + 1) * (J + 2)))
  c.u <- a.u * (i - (J + 1) / 2)^2 + b.u
  a.u.ref <- 6 * sqrt(5 / (J.ref^5 - 5 * J.ref^3 + 4 * J.ref))
  b.u.ref <- -1 / 2 * (sqrt(5) * J.ref^2 - sqrt(5)) / 
       (sqrt((J.ref - 2) * (J.ref - 1) * J.ref * (J.ref + 1) * (J.ref + 2)))
  c.u.ref <- a.u.ref * (i.ref - (J.ref + 1) / 2)^2 + b.u.ref

  # calculate rank histograms
  rh.ens <- rankhist(ens, ver)
  rh.ref <- rankhist(ens.ref, ver)

  # calculate derived quantities for chi^2 tests
  xi <- (rh.ens - e.i) / sqrt(e.i)
  xi.ref <- (rh.ref - e.i.ref) / sqrt(e.i.ref)
  
  # calculate the score differences
  score.diffs <- c(
    # difference in pearson chi^2 statistic
    sum(xi.ref * xi.ref) - sum(xi * xi),
    # difference in jolliffe-primo slope statistic
    sum(xi.ref * c.lin.ref)^2 - sum(xi * c.lin)^2,
    # difference in jolliffe-primo Ushape statistic
    sum(xi.ref * c.u.ref)^2 - sum(xi * c.u)^2
  )

  # bootstrap the null distribution of "score" differences
  ens.combi <- cbind(ens, ens.ref)
  s.H0 <- t(replicate(n.boot, {
    ens.shuf <- ens.combi[, sample(1:(K+K.ref), K+K.ref, replace=TRUE)]
    ens.shuf <- ens.shuf[sample(1:N, N, replace=TRUE), ]
    rh1 <- rankhist(ens.shuf[,1:K], ver)
    x1 <- (rh1 - e.i) / sqrt(e.i)
    rh2 <- rankhist(ens.shuf[,(K+1):(K+K.ref)], ver)
    x2 <- (rh2 - e.i.ref) / sqrt(e.i.ref)
    c(sum(x2 * x2) - sum(x1 * x1),
      sum(x2 * c.lin.ref)^2 - sum(x1 * c.lin)^2,
      sum(x2 * c.u.ref)^2 - sum(x1 * c.u)^2)
  }))

  # calculate p values of observed statistics under the null
  p.values <- c(
    1 - ecdf(s.H0[, 1])(score.diffs[1]),
    1 - ecdf(s.H0[, 2])(score.diffs[2]),
    1 - ecdf(s.H0[, 3])(score.diffs[3])
  )
  
  # calculate ranks of ens and ens.ver
  ranks.ens <- apply(cbind(ver, ens), 1, rank, ties.method="random")[1, ]
  ranks.ref <- apply(cbind(ver, ens.ref), 1, rank, ties.method="random")[1, ]
  r.df <- data.frame(ranks.ens=ranks.ens, ranks.ref=ranks.ref)

  # bootstrap the sampling distribution of the "score" differences
  s.H1 <- t(replicate(n.boot, {
    ens.shuf <- ens[, sample(1:K, K, replace=TRUE)]
    ens.ref.shuf <- ens.ref[, sample(1:K.ref, K.ref, replace=TRUE)]
    inds <- sample(1:N, N, replace=TRUE)
    rh.ens <- rankhist(ens.shuf[inds, ], ver[inds])
    rh.ref <- rankhist(ens.ref.shuf[inds, ], ver[inds])
    # calculate necessary quantities
    xi <- (rh.ens - e.i) / sqrt(e.i)
    xi.ref <- (rh.ref - e.i.ref) / sqrt(e.i.ref)
    # difference in pearson chi^2 statistic
    prsn.diff <- sum(xi.ref * xi.ref) - sum(xi * xi) 
    # difference in jolliffe-primo linear
    jplin.diff <- sum(xi.ref * c.lin.ref)^2 - sum(xi * c.lin)^2
    # difference in jolliffe-primo U
    jpu.diff <- sum(xi.ref * c.u.ref)^2 - sum(xi * c.u)^2
    # return
    c(prsn.diff, jplin.diff, jpu.diff)
  }))

  # estimate some sampling quantiles of the score differences
  quantls <- rbind(
    quantile(s.H1[, 1], probs=c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99)),
    quantile(s.H1[, 2], probs=c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99)),
    quantile(s.H1[, 3], probs=c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99))
  )

  # return
  ret.df <- cbind(score.diffs, p.values, quantls)
  colnames(ret.df) <- c("score.diff", "p.value", paste0("Q", c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99)))
  rownames(ret.df) <- c("pearson.chi2", "jp.slope", "jp.convex")
  ret.df
}


################################
#
# MINIMUM SPANNING TREE RANK HISTOGRAM
#
#  *** WORK IN PROGRESS, USE AT YOUR OWN RISK ***
#
# ens ... array of dim N, K, D
#         N is the number of forecasts, there are K ensemble members
#         per forecast, each ensemble member is a vector of dimension D;
#         two ensembles (N*K matrices) can be pieced together along a new
#         dimension D by abind(ens1, ens2, along=3) to yield a multidimensional
#         ensemble
# ver ... array of dim N, D; the corresponding D-dimensional
#         verifications
#
# references: Smith & Hansen (2007)
#             Wilks (2007)
#
################################
mstrankhist <- function(ens, ver) {
  # ens[i,,] is an K * D matrix, the rows are the ensemble members
  N <- dim(ens)[1]
  K <- dim(ens)[2]
  D <- dim(ens)[3]
  r <- sapply(1:N, function(i) {
         e <- rbind(ens[i, , ], ver[i, ])
         K <- nrow(e)
         v <- sapply(seq(K), function(i) {
                x <- e[-i, ]
                d <- dist(x, diag=TRUE)
                sum(as.matrix(d) * mst(d)) / 2.0
              })
         tail(rank(v, ties.method="random"), 1)
       })
  return(r)
}

################################
#                              #
# PLOT RANK HISTOGRAM          #
#                              #
################################
PlotRankhist <- function(rank.hist, mode=c("raw", "prob.paper")) {
#
#
# Plot the rank histogram raw or on probability paper
#
# Usage: PlotRankhist(rh, mode="raw")
#
# Arguments:
#   * rank.hist ... a vector of rank counts (see function `rankhist()`)
#   * mode      ... whether to plot raw or on probability paper to 
#                   assess flatness
#
# Details:
#
#   * the `raw` mode simply plots a barplot of the rank histogram counts
#   * the `prob.paper` mode transforms the observed rank histogram counts to
#     cumulative probabilities under Binomial(N, 1/(K+1)), plots them on a logit
#     scale, and adds simultaneous consistency intervals
#
# Return value:
#   * none (yet)
#
# Author: 
#
#    Stefan Siegert 
#    s.siegert@exeter.ac.uk 
#    October 2013
#
# Example:
# 
#   ens <- matrix(rnorm(500), 100, 5)
#   ver <- rnorm(100)
#   rh <- rankhist(ens, ver)
#   PlotRankhist(rank.hist = rh, mode="prob.paper")
#
# References: 
#   Broecker 2008
#
################################
  if (mode == "prob.paper") {
    N <- sum(rank.hist)
    K <- length(rank.hist) - 1
    # cumulative binomial likelihood of the observed rank counts
    nuh <- pbinom(q=rank.hist, size=N, prob=1/(K+1))
    # log-odds ratio
    lornuh <- log(nuh / (1 - nuh))
    # clip very small and large bars
    clipval <- 8
    i.clip.min <- which(lornuh < -clipval)
    lornuh[i.clip.min] <- -clipval
    i.clip.max <- which(lornuh >  clipval)
    lornuh[i.clip.max] <-  clipval
    # prepare for plotting
    offs <- 8
    bar.wd <- 0.9
    par(oma=c(0, 0, 0, 3), cex.lab=0.8, cex.axis=0.8)
    plot(NULL, xlim=c(0,K+2), ylim=c(0,2*offs), axes=F, xlab="rank i", ylab=expression(nu[i]))
    b <- barplot(lornuh + offs, add=T, axes=FALSE, width=bar.wd, col=gray(0.5))
    points(b[i.clip.min], lornuh[i.clip.min]+offs, pch=25, bg="black", cex=.6)
    points(b[i.clip.max], lornuh[i.clip.max]+offs, pch=24, bg="black", cex=.6)
    # axis labels
    xlabels <- matrix(t(cbind(paste(seq(1,K+1,2)),"")),nrow=1)[1:(K+1)]
    axis(1, at=b, labels=xlabels)
    pvals <- c(0.001,0.01,0.1,0.5,0.9,0.99,0.999)
    yvals <- log(pvals/(1-pvals))
    axis(2, at=yvals+offs, labels=pvals, las=2)
    # confidence intervals, corrected for multiple testing
    ci.str <- c("90%","95%","99%")
    k <- 0
    for (p in c(0.9,0.95,0.99)) {
      k <- k+1
      del.b <- diff(b)[1]
      ci <- c(1-p^(1/(K+1)),p^(1/(K+1)))
      ci <- log(ci/(1-ci))+offs
      axis(4,at=ci,labels=c("",""),tcl=0.5,line=k)
      axis(4,at=log(c(0.1,0.9)/(1-c(0.1,0.9)))+offs,labels=c("",""),col="white",lwd=3,line=k)
      mtext(ci.str[k],side=4,padj=-1,line=k,cex=0.8)
      lines(x=c(b[1]-.7*bar.wd,b[K+1]+0.7*bar.wd),y=rep(ci[1],2),lty=2,xpd=TRUE)
      lines(x=c(b[1]-.7*bar.wd,b[K+1]+0.7*bar.wd),y=rep(ci[2],2),lty=2,xpd=TRUE)
    }
  } else if (mode == "raw") {
    par(mar=c(3, 3, 1, 0), cex.lab=0.8, cex.axis=0.8)
    bp <- barplot(rank.hist, xlab=NA, ylab=NA, col=gray(0.5))
    axis(1, at=bp, line=.2, labels=paste(1:length(rank.hist)))
    mtext(side=1, text="rank i", line=2, cex=.8)
    mtext(side=2, text="count", line=2, cex=.8)
  }
}



