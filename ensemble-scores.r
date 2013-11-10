################################
#
# THE FAIR BRIER SCORE
#
# ens ... ensemble matrix (N*K)
# ver ... verification vector (N)
# tau ... exceedance threshold that defines the event
#
# reference: Ferro (2013) "Fair scores for ensemble forecasts"
#
################################
fairbrier <- function(ens, ver, tau=0.5) {
  ver <- matrix(ver, ncol=1)
  if (is.null(dim(ens))) {
    ens <- matrix(ens, nrow=1)
  }
  stopifnot(nrow(ens) == nrow(ver))

  K <- ncol(ens)
  i <- rowSums(ens > tau)
  j <- 1 * (ver > tau)
  fb <- (j - i / K) ^ 2 - i * (K - i) / K / K / (K - 1)
  return(fb)
}

################################
#
# ANALYZE DIFFERENCE IN THE FAIR BRIER SCORE BETWEEN TWO ENSEMBLE
# FORECASTING SYSTEMS FOR THE SAME VERIFICATION
#
# ens     ... ensemble to be tested
# ens.ref ... reference forecast ensemble
# ver     ... verifications
# tau     ... threshold, whose exceedance defines the event
# n.boot  ... number of bootstrap samples
#
################################
AnalyzeFairBrierDifference <- 
function(ens, ens.ref, ver, tau=0.5, n.boot=500) {

  ver <- matrix(ver, ncol=1)
  if (is.null(dim(ens))) {
    ens <- matrix(ens, nrow=1)
    ens.ref <- matrix(ens.ref, nrow=1)
  }
  stopifnot(all(c(nrow(ens), nrow(ens.ref)) == nrow(ver)))

  K <- ncol(ens)
  K.ref <- ncol(ens.ref)

  # calculate fair Brier score difference
  br.ens <- fairbrier(ens, ver, tau)
  br.ref <- fairbrier(ens.ref, ver, tau)
  br.diff <- mean(br.ref - br.ens)

  # bootstrap the null distribution and estimate p-value
  ens.combi <- cbind(ens, ens.ref)
  s.H0 <- replicate(n.boot, {
    ens.shuf <- ens.combi[, sample(1:(K+K.ref), K+K.ref)]
    br <- fairbrier(ens.shuf[,1:K], ver, tau)
    br.ref <- fairbrier(ens.shuf[,(K+1):(K+K.ref)], ver, tau)
    mean(br.ref - br)
  })
  p.value <- 1 - ecdf(s.H0)(br.diff)

  # bootstrap the sampling distribution and estimate some quantiles
  br.df <- data.frame(br.ens=br.ens, br.ref=br.ref)
  f <- function(br.df, inds) {
    with(br.df, mean(br.ref[inds] - br.ens[inds]))
  }
  test.out <- boot(data=br.df, statistic=f, R=n.boot)
  quantls <- quantile(test.out[["t"]], c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99))
  #return
  ret.df <- c(br.diff, p.value, quantls)
  names(ret.df) <- c("fair.brier.diff", "p.value", paste0("Q", c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99)))
  return(ret.df)
}


################################
#
# THE FAIR CONTINUOUSLY RANKED PROBABILITY SCORE
#
# ens ... ensemble values 
# ver ... verification
#
# references: * Gneiting, Raftery (2007) "Probabilistic forecasts,
#               calibration and sharpness"
#             * Ferro, Richardson, Weigel (2008) "On the effect 
#               of ensemble size on the discrete and continuous
#               ranked probability scores"
#
################################
faircrps <- 
function(ens, ver) {
  ver <- as.vector(ver)
  if (length(ver) == 1) {
  # single instance
    if (is.matrix(ens)) {
      ens <- as.vector(ens)
    }
    K <- length(ens)
    if (K == 1) {
    # for one ensemble member, the crps reduces to the absolute error
      crps <- abs(ens - ver)
    } else {
      crps <- mean(abs(ens - ver)) - sum(dist(ens)) / K / (K - 1)
    }
  } else {
  # multiple instances
    N <- length(ver)
    K <- ncol(ens)
    stopifnot(length(ver) == nrow(ens))
    if (K == 1) {
      crps <- abs(ens - ver)
    } else {
      crps <- sapply(1:N, function(i) 
                mean(abs(ens[i,] - ver[i])) - sum(dist(ens[i,])) / K / (K - 1)
              )
    }
  }
  return(crps)
}
  


################################
#
# ANALYZE DIFFERENCE IN THE FAIR CRPS BETWEEN TWO ENSEMBLE
# FORECASTING SYSTEMS FOR THE SAME VERIFICATION
#
# ens     ... number of ensemble members that predict the 
#              event in the ensemble
# ens.ref ... number of ensemble members that predict the 
#            event in the reference ensemble
# ver     ... equals one if the event happens, zero otherwise
# test    ... t.test, bootstrapping
# size    ... size of the test, probability of false detection
#
################################
AnalyzeFairCrpsDifference <- 
function(ens, ens.ref, ver, n.boot=500) {

  ver <- matrix(ver, ncol=1)
  if (is.null(dim(ens))) {
    ens <- matrix(ens, nrow=1)
    ens.ref <- matrix(ens.ref, nrow=1)
  }
  stopifnot(all(c(nrow(ens), nrow(ens.ref)) == nrow(ver)))

  K <- ncol(ens)
  K.ref <- ncol(ens.ref)

  # calculate fair crps difference
  crps.ens <- faircrps(ens, ver)
  crps.ref <- faircrps(ens.ref, ver)
  crps.diff <- mean(crps.ref - crps.ens)

  # bootstrap the null distribution and estimate p-value
  ens.combi <- cbind(ens, ens.ref)
  s.H0 <- replicate(n.boot, {
    ens.shuf <- ens.combi[, sample(1:(K+K.ref), K+K.ref)]
    crps <- faircrps(ens.shuf[,1:K], ver)
    crps.ref <- faircrps(ens.shuf[,(K+1):(K+K.ref)], ver)
    mean(crps.ref - crps)
  })
  p.value <- 1 - ecdf(s.H0)(crps.diff)

  # bootstrap the sampling distribution and estimate some quantiles
  crps.df <- data.frame(crps.ens=crps.ens, crps.ref=crps.ref)
  f <- function(crps.df, inds) {
    with(crps.df, mean(crps.ref[inds] - crps.ens[inds]))
  }
  test.out <- boot(data=crps.df, statistic=f, R=n.boot)
  quantls <- quantile(test.out[["t"]], c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99))
  #return
  ret.df <- c(crps.diff, p.value, quantls)
  names(ret.df) <- c("fair.crps.diff", "p.value", paste0("Q", c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99)))
  return(ret.df)
}
