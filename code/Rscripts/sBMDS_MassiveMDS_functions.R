library(MassiveMDS)

##### target functions
# (the log-posterior which is proportional to log-likelihood + log-priors)

s_target_no_sigma_prior <- function(locations, engine) {
  # theta is the latent variable
  output <- MassiveMDS::getLogLikelihood(engine) +
    # independent, Gaussian prior for theta centered at 0 & sd = 1
    sum(mvtnorm::dmvnorm(locations, mean = rep(0, engine$embeddingDimension),
                         sigma = diag(rep(1, engine$embeddingDimension)), log = TRUE))
  return(output)
}

s_target <- function(locations, engine) {
  sigmasq <- 1 / engine$precision
  output <- MassiveMDS::getLogLikelihood(engine) +
    # independent, Gaussian prior for theta centered at 0 & sd = 1
    sum(mvtnorm::dmvnorm(locations, mean = rep(0, engine$embeddingDimension),
                         sigma = diag(rep(1, engine$embeddingDimension)), log = TRUE)) +
    # inverse-gamma prior for sigma^2
    dinvgamma(sigmasq, shape = 1, rate = 1, log = TRUE)
  return(output)
}

### gradient for target function

sbmds_grad <- function(locations, engine){
  gradient <- MassiveMDS::getGradient(engine) - locations
  return(gradient)
}

##### delta function, rate of sd stepsize increase for adaptive MH/HMC

delta <- function(n) {
  return( min(0.01, n^(-0.5)) )
}

##### non-adaptive sigma, adaptive latent variable, Metropolis

sbmds_metropolis <- function(maxIts, data, targetAccept = 0.8, stepSize = 1) {

  sigmasq <- 1 / engine$precision
  dims <- engine$embeddingDimension
  n <- engine$locationCount

  # create the chain
  chain <- array(0, dim = c(maxIts, n, dims))

  # specify the first random value
  chain[1, , ] <- mvtnorm::rmvnorm(n, mean = rep(0, dims), sigma = diag(dims))

  totalAccept <- rep(0, maxIts)
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 50   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)

  for (s in 2:maxIts) {

    thetaStar <- matrix(rnorm(n * dims, mean = chain[s - 1, , ], sd = stepSize),
                        ncol = dims, nrow = n)
    # smaller proposal variance = more acceptances
    # proposal: normal dist
    u <- runif(1)
    # Metropolis, A = target(thetaStar)/target(previous iteration)
    # target on log scale
    logA <- s_target_no_sigma_prior(thetaStar, engine) -
      s_target_no_sigma_prior(chain[s - 1, , ], engine)

    if(log(u) < logA) {
      chain[s, , ] <- thetaStar # ACCEPT !! # next iteration to thetastar
      totalAccept[s] <- 1
      Acceptances = Acceptances + 1 # acceptance counter
    } else {
      chain[s, , ] <- chain[s - 1, , ] # REJECT !!
      # next iteration stays at same iteration
    }

    SampCount <- SampCount + 1

    # tune
    if (SampCount == SampBound) {
      AcceptRatio <- Acceptances / SampBound
      if (AcceptRatio > targetAccept) {
        stepSize <- stepSize * (1 + delta(s - 1)) # increase stepsize
      } else {
        stepSize <- stepSize * (1 - delta(s - 1)) # decrease stepsize
      }

      # reset Sampcount and Acceptances
      SampCount <- 0
      Acceptances <- 0
    }

    if (s %% 100 == 0) cat("Iteration ", s, "\n","stepSize: ", stepSize, "\n")
  }

  cat("Acceptance rate: ", sum(totalAccept)/(maxIts - 1))

  return(chain)
}



