####################################################
#                                                  #
#  FAIR RELIABILITY DIAGRAM FOR ENSEMBLE FORECASTS #
#                                                  #
####################################################
fair.rel.diag <- function(i, j, K, nboot=1000, bar.probs=c(0.025, 0.975), 
                     plot=FALSE, plot.refin=TRUE) 
  # Plot reliability diagram for an ensemble forecast
  #
  # Usage: fair.rel.diag(i, j, K, nboot, bar.probs, plot, plot.refin)
  #
  # Arguments:
  #
  #    i ... vector of length N, i[k] has the number of ensemble members that predict
  #          the event j[k], k = 1...N
  #    j ... j[k] = 1 if the event happened at instance k, j[k] = 0 otherwise
  #    K ... number of ensemble members, i[k] element {0...K} for all k = 1...N
  #    nboot ... number of bootstrap resamples for estimating consistency bars; 
  #              if nboot==0, NAs are returned as consistency bars
  #    bar.probs ... vector of length 2, lower and upper confidence limit
  #    plot ... boolean; whether to plot the reliability diagram
  #    plot.refin ... boolean; whether to include the refinement diagram in the lower right corner
  #
  # Return value:
  #
  #    a data frame of K+1 rows with the following columns:
  #
  #      * i          ... (0, 1, ..., K)
  #      * cond.probs ... observed conditional frequency of event, given i
  #      * H0.line    ... for a reliable ensemble, the cond.probs[i] = H0.line[i] for all i
  #      * cbar.lo    ... lower limit consistency of consistency bar[i], as specified by user
  #      * cbar.hi    ... upper limit consistency of consistency bar[i], as specified by user
  #
  # Author: 
  #
  #    Stefan Siegert 
  #    s.siegert@exeter.ac.uk 
  #    October 2013
  #
  # Example:
  #
  #    N <- 1000
  #    K <- 10
  #    mu <- runif(N)
  #    ens <- matrix(rnorm(N*K, mean=mu, sd=1), N, K)
  #    ver <- rnorm(N, mean=mu, sd=1)
  #    tau <- 1
  #    i <- rowSums(ens > tau)
  #    j <- 1 * (ver > tau)
  #    rd <- fair.rel.diag(i, j, K, plot=TRUE)
  #    print(rd)
{
  # parameters and variables
  N <- length(j)
  p.ens <- i / K

  # estimate hyper-parameters (method of moments)
  m1 <- mean(i)
  m2 <- mean(i*i)
  a <- (K * m1 - m2) / (K * (m2 / m1 - m1 - 1) + m1)
  b <- (K - m1) * (K - m2 / m1) / (K * (m2 / m1 - m1 - 1) + m1)

  # if either parameter is negative, use maximum-likelihood estimation
  if (a < 0 | b < 0) {
    nloglik <- function(par,i,K) {
      a <- par[1]
      b <- par[2]
      # return negative log likelihood
      # terms that dont depend on a or b were removed
      -sum(lbeta(i+a, K-i+b) - lbeta(a,b))
    }
    loglik.opt <- optim(par=c(1,1), fn=nloglik, i=i, K=K)
    stopifnot(loglik.opt[["convergence"]] == 0)
    a <- loglik.opt[["par"]][1]
    b <- loglik.opt[["par"]][2]
  }

  # estimate calibration function
  nbins <- K + 1
  brx <- seq(0, 1, length.out=nbins+1) +
         c(-.Machine$double.eps, rep(0, nbins-1), .Machine$double.eps)
  h <- hist(p.ens, breaks=brx, plot=FALSE)$counts
  g <- hist(p.ens[j==1], breaks=brx, plot=FALSE)$counts
  obar.i <- g / h

  # consistency bars by resampling
  if (nboot) {
    resamp.mat <- matrix(nrow=0, ncol=nbins)
    for (ii in 1:nboot) {
      p.hat <- rbeta(N, a, b)
      y.hat <- rbinom(N, 1, p.hat)
      k.hat <- rbinom(N, K, p.hat)
      p.ens.hat <- k.hat / K
      h <- hist(p.ens.hat, breaks=brx, plot=FALSE)$counts
      g <- hist(p.ens.hat[y.hat==1], breaks=brx, plot=FALSE)$counts
      resamp.mat <- rbind(resamp.mat, g / h)
    }
    cons.bars <- apply(resamp.mat, 2, quantile, probs=bar.probs, na.rm=TRUE)
  } else {
    cons.bars <- matrix(NA, ncol=K+1, nrow=2)
  }

  # return object
  ret.df <- data.frame(i=(0:K), cond.probs=obar.i, H0.line=(a+(0:K))/(a+b+K), 
                       cbar.lo=cons.bars[1,],cbar.hi=cons.bars[2,])

  # plot
  if (plot) {
    old.par <- par(no.readonly = TRUE) 
    on.exit(par(old.par))
    plot(NULL, xlim=c(0,K), ylim=c(0,1), xlab="i", ylab="observed relative frequencies")
    with(ret.df, {
      for (ii in i) lines(rep(ii,2), c(cbar.lo[ii+1], cbar.hi[ii+1]), col=gray(0.7), lwd=5)
      lines(i, H0.line, lty=2)
      points(i, cond.probs, pch=1, lwd=2)
    })
    if (plot.refin) {
      # refinement histogram and beta-binomial fit in lower corner
      k.hist <- hist(i, breaks=seq(0,K+1)-0.5, plot=FALSE)$counts
      k.fit <- N * beta((0:K)+a, K-(0:K)+b) * choose(K, (0:K)) / beta(a, b)
      pp<- par("plt")
      par("plt" = c(pp[2] - 0.2 , pp[2],  pp[3], pp[3]+ 0.2) )
      par(new = TRUE)
      bpl <- barplot(k.hist, axes=FALSE, ylim=c(0, max(k.hist)*1.2), border=NA)
      points(bpl, k.fit, pch=15, cex=0.6, type="b")
      axis(4, at=pretty(c(0,max(k.hist)), 2), las=2, cex.axis=.8)
      box()
    }
  }

  return(ret.df)
}




######################################################################
#                                                                    #
# RELIABILITY DIAGRAM FOR A COLLECTION OF PROBABILITY FORECASTS      #
#                                                                    #
######################################################################
rel.diag <- function(probs, ver, nbins=10, nboot=1000, plot=FALSE, plot.refin=TRUE) {
  #
  # Plot reliability diagram for a probability forecast
  #
  # Usage: rel.diag(probs, ver, nbins, nboot)
  #
  # Arguments:
  #
  #    probs ... vector of length N, probs[k] has the predicted probability for the event ver[k] 
  #    ver ... ver[k] = 1 if the event happened at instance k, ver[k] = 0 otherwise
  #    nbins ... number of bins to discretize the forecast probabilities
  #    nboot ... number of bootstrap resamples for estimating consistency bars
  #              if nboot==0, NAs are returned as consistency bars
  #    plot ... boolean; whether to plot the reliability diagram
  #    plot.refin ... boolean; whether to plot the small refinement histogram in lower right corner
  #
  # Return value:
  #
  #    a data frame of K+1 rows with the following columns:
  #
  #      * p.avgs     ... in-bin averages of the forecast probabilities
  #      * cond.probs ... observed conditional frequency of event, given i
  #      * cbar.lo    ... lower limit consistency of consistency bar[i], as specified by user
  #      * cbar.hi    ... upper limit consistency of consistency bar[i], as specified by user
  #
  # Author: 
  #
  #    Stefan Siegert 
  #    s.siegert@exeter.ac.uk 
  #    October 2013
  #
  # Example:
  #
  #    N <- 1000
  #    p <- rbeta(N, 1, 3)
  #    y <- rbinom(N, 1, p)
  #    rd <- rel.diag(p, y, plot=TRUE)
  #    print(rd)
  #
  #
  # change log:
  #
  #
  #  2013/10/31:
  #  * return summary data as data frame
  #  * added options `plot` and `plot.refin`
  #
  #  2013/08/20:
  #  * points are plotted at in-bin-averages, not at bin centres
  #  * legend has been removed
  #  * consistency bars have been added, calculated by a resampling technique
  #  * see Broecker (2007) http://dx.doi.org/10.1175/WAF993.1 for details
  #  * the bars are pointwise 2.5% ... 97.5% intervals around the hypothesis of reliability
  #  * dependency on package "verification" was removed
  #
  # Author: Stefan Siegert <s.siegert@exeter.ac.uk>
  #
  # based on previous version by Caio Coelho and the routine 
  # reliability.plot.default of the R-package `verification`
  #


  n <- length(ver)

  # estimate refinement function
  brx <- seq(0, 1, length.out=nbins+1) + 
         c(-.Machine$double.eps, rep(0, nbins-1), .Machine$double.eps)
  h <- hist(probs, breaks=brx, plot=FALSE)$counts        

  # estimate calibration function
  g <- hist(probs[ver==1], breaks=brx, plot=FALSE)$counts
  obar.i <- g / h 
  
  # calculate in-bin averages
  p.bins <- as.numeric(cut(probs, breaks=brx))
  p.avgs <- sapply(seq(nbins), 
                   function(ii) mean(probs[p.bins == ii], na.rm=TRUE))

  # consistency resampling (broecker and smith 2007)
  resamp.mat <- matrix(nrow=0, ncol=nbins)
  if (nboot) {
    for (i in 1:nboot) {
      p.hat <- sample(x=probs, size=n, replace=TRUE)
      x.hat <- rbinom(n=n, size=1, prob=p.hat)
      hh <- hist(p.hat, breaks=brx, plot=FALSE)$counts        
      gg <- hist(p.hat[x.hat==1], breaks=brx, plot=FALSE)$counts
      resamp.mat <- rbind(resamp.mat, gg / hh)
    }
    cons.bars <- apply(resamp.mat, 2, 
                       function(z) quantile(z, c(.025, .975), na.rm=TRUE))
  } else {
    cons.bars <- matrix(NA, ncol=nbins, nrow=2)
  }

  if (plot) {
    # reliability plot
    old.par <- par(no.readonly = TRUE) 
    on.exit(par(old.par))
    plot(NULL, xlim = c(0,1), ylim = c(0,1),
       xlab= "Forecast probability",
       ylab="Observed relative frequency")
    # consistency bars
    for (i in 1:length(p.avgs)) {
        lines(c(p.avgs[i], p.avgs[i]), cons.bars[, i], col="#CCCCCC", lwd=6)
    }
    # reliability points and diagonal
    points(p.avgs, obar.i, col = "black", pch = 1, lwd=2)
    abline(0,1)
    if (plot.refin) {
      # refinement histogram in lower corner
      pp<- par("plt")
      par("plt" = c(pp[2] - 0.2 , pp[2],  pp[3], pp[3]+ 0.2) )
      par(new = TRUE)
      barplot(h, axes = FALSE, axisnames = FALSE)
      #hist(probs, axes=FALSE, main=NA, xlab=NA, ylab=NA, breaks=brx, col="gray")
      #lines(density(probs, width=0.01), col="black")
      axis(4)
      box() 
    }
  }

  # return object
  ret.df <- data.frame(p.avgs=p.avgs, cond.probs=obar.i, 
                       cbar.lo=cons.bars[1,], cbar.hi=cons.bars[2,])
  return(ret.df)
}

