stopifnot(require(boot))  # bootstrapping
#stopifnot(require(ape))   # minimum spanning tree
#stopifnot(require(abind)) # multidimensional arrays 

###########################
#                         #
# CREATE TOY ENSEMBLE AND #
#    VERIFICATION         #
#                         #
###########################
GenerateToyData <- function(mu.ens=0, sd.ens=1, 
                            mu.ver=0, sd.ver=1, 
                            mu.ref=NA, sd.ref=NA,
                            K=10, K.ref=10, N=100) 
#
# Generate artificial Gaussian ensemble data and verifications
#
# Usage: 
#   td <- GenerateToyData()
#
# Arguments:
#
#   N ... number of instances
#   K ... number of ensemble members of `ens`
#   K.ref ... number of ensemble members of `ens.ref`
#   mu.ver, sd.ver ... means and variances of the verification
#   mu.ens, sd.ens ... means and variances of the ensemble
#   mu.ref, sd.ref ... means and variances of the reference ensemble
#
# Return value:
#   * a list containing two or three objects, depending on whether any of `mu.ref` or `sd.ref` is NA:
#     + ver ... a vector of Gaussian variables drawn from N(mu.ver, sd.ver)
#     + ens ... a N*K matrix, each row is an ensemble drawn from N(mu.ens, sd.ens)
#     + ens.ref ...  a N*K.ref matrix, each row is an ensemble drawn from
#       N(mu.ref, sd.ref); only returned if neither of mu.ref and sd.ref is NA
#
# Author:
#   Stefan Siegert 
#   s.siegert@exeter.ac.uk 
#   October 2013
#
# Example:
#   mu.ver <- runif(100)
#   mu.ens <- mu.ver + runif(100, -0.1, .1)
#   mu.ref <- mu.ver + runif(100, -0.2, .2)
#   td <- GenerateToyData(mu.ver=mu.ver, mu.ens=mu.ens, mu.ref=mu.ref,
#                         sd.ref=1, K=10, K.ref=10)
#
#
{
  if (is.na(mu.ref) | is.na(sd.ref)) {
    l <- list(
        ens = matrix(rnorm(N*K, mu.ens, sd.ens), N, K),
        ver = rnorm(N, mu.ver, sd.ver))
  } else {
    l <- list(
        ens = matrix(rnorm(N*K, mu.ens, sd.ens), N, K),
        ens.ref = matrix(rnorm(N*K.ref, mu.ref, sd.ref), N, K.ref),
        ver = rnorm(N, mu.ver, sd.ver))
  }
  return(l)
}



