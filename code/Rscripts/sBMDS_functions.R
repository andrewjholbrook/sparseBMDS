##### raw likelihood

## shows difference in computational time as a function of band.no
sumr <- function(latent, D, sigmasq, band.no){
  D_latent <- as.matrix(dist(latent, diag = TRUE, upper = TRUE))
  n <- dim(D)[1]
  SSR_val = 0
  logsum = 0
  pad <- matrix(NA, nrow = n, ncol = band.no)
  D_latent_pad <- cbind(D_latent, pad)
  D_pad <- cbind(D, pad)

  for (i in 1:n){
    for (j in 1:band.no){
      SSR_val = sum(SSR_val, (D_pad[i, i + j] - D_latent_pad[i, i + j])^2,
                    na.rm = TRUE)
      stat = D_latent_pad[i, i + j] / sqrt(sigmasq)
      logsum = sum(logsum, log(pnorm(stat)), na.rm = TRUE)
    }
  }
  r = SSR_val/(2 * sigmasq) + logsum
  return(r)
}

sbmds_ll <- function(latent, D, sigmasq, band.no){
  # latent variables (nxk), observed distance matrix D (nxn), sigma^2
  # D_latent <- as.matrix(dist(latent, diag = TRUE, upper = TRUE))
  n <- dim(D)[1]
  m <- n * band.no - band.no * (band.no + 1) / 2
  term1 <- (m / 2) * log(sigmasq)
  term2 <- sumr(latent, D, sigmasq, band.no)
  loglik <- -(term1 + term2)
  return(loglik)
}

## fast version
sumr_fast <- function(latent, D, sigmasq, band.no){
  D_latent <- as.matrix(dist(latent, diag = TRUE, upper = TRUE))
  n <- dim(D)[1]
  band_latent <- matrix(NA, nrow = n - 1, ncol = band.no)
  band_obs <- matrix(NA, nrow = n - 1, ncol = band.no)
  pad <- matrix(NA, nrow = n, ncol = band.no)
  D_latent_pad <- cbind(D_latent, pad)
  D_pad <- cbind(D, pad)

  for (i in 1:(n - 1)){
    band_latent[i, ] <- D_latent_pad[i, i + 1:band.no]
    band_obs[i, ] <- D_pad[i, i + 1:band.no]
  }

  SSR_val = (band_obs - band_latent)^2
  stat = band_latent / sqrt(sigmasq)
  logsum = log(pnorm(stat))
  r = sum(SSR_val/(2 * sigmasq), logsum, na.rm = TRUE)
  return(r)
}

sbmds_ll_fast <- function(latent, D, sigmasq, band.no){
  # latent variables (nxk), observed distance matrix D (nxn), sigma^2
  # D_latent <- as.matrix(dist(latent, diag = TRUE, upper = TRUE))
  n <- dim(D)[1]
  m <- n * band.no - band.no * (band.no + 1) / 2
  term1 <- (m / 2) * log(sigmasq)
  term2 <- sumr_fast(latent, D, sigmasq, band.no)
  loglik <- -(term1 + term2)
  return(loglik)
}

##### raw gradient

sbmds_grad <- function(latent, D, sigmasq, band.no) {
  # padding to prevent errors
  dims <- dim(latent)[2]
  n <- dim(latent)[1]
  D_latent <- as.matrix(dist(latent, diag = TRUE, upper = TRUE))

  pad <- matrix(NA, nrow = n, ncol = band.no) # add padding column-wise
  pad2 <- matrix(NA, nrow = band.no, ncol = dims) # add padding row-wise

  # padding to the left & right of dist matrices
  D_latent_pad <- cbind(pad, D_latent, pad)
  D_pad <- cbind(pad, D, pad)
  # padding above & below latent variables
  latent_pad <- rbind(pad2, latent, pad2)
  # computing pnorm & dnorm to entire D_latent matrix
  D_latent_pnorm <- pnorm(D_latent_pad / sqrt(sigmasq))
  D_latent_dnorm <- dnorm(D_latent_pad / sqrt(sigmasq))

  grad_sll <- matrix(NA, nrow = n, ncol = dims)

  for (i in 1:n){
    step3 <- matrix(0, nrow = 2 * band.no, ncol = 2)
    for (j in 1:band.no){
      # accounting for couplings to the left & right of d_ii
      step1 <- (D_latent_pad[i, c(i + band.no - j, i + band.no + j)] -
                    D_pad[i, c(i + band.no - j, i + band.no + j)]) / sigmasq
      step2 <- D_latent_dnorm[i, c(i + band.no - j, i + band.no + j)] /
        (sqrt(sigmasq) * D_latent_pnorm[i, c(i + band.no - j, i + band.no + j)])
      # elements to the left
      step3[j, ] <- (step1 + step2)[1] *
        ((latent[i, ] - latent_pad[i + band.no - j, ]) /
           D_latent_pad[i, i + band.no - j])
      # elements to the right
      step3[j + band.no, ] <- (step1 + step2)[2] *
        ((latent[i, ] - latent_pad[i + band.no + j, ]) /
           D_latent_pad[i, i + band.no + j])
    }

    # combine both types of couplings
    grad_sll[i, ] <- - colSums(step3, na.rm = TRUE) - latent[i, ]
  }
  return(grad_sll)
}

##### target functions
# (the log-posterior which is proportional to log-likelihood + log-priors)

s_target_no_sigma_prior <- function(theta, D, sigmasq, band.no) {
  dims <- dim(theta)[2]
  # theta is the latent variable
  output <- sbmds_ll(theta, D, sigmasq, band.no) +
    # independent, Gaussian prior for theta centered at 0 & large sd
    sum(mvtnorm::dmvnorm(theta, mean = rep(0, dims),
                         sigma = diag(rep(1, dims)), log = TRUE))
  return(output)
}

s_target_no_sigma_prior_fast <- function(theta, D, sigmasq, band.no) {
  dims <- dim(theta)[2]
  output <- sbmds_ll_fast(theta, D, sigmasq, band.no) +
    sum(mvtnorm::dmvnorm(theta, mean = rep(0, dims),
                         sigma = diag(rep(1, dims)), log = TRUE))
  return(output)
}

s_target <- function(theta, D, sigmasq, band.no) {
  output <- sbmds_ll(theta, D, sigmasq, band.no) +
    # independent, standard Gaussian prior for theta
    sum(mvtnorm::dmvnorm(theta, log = TRUE)) +
    # inverse-gamma prior for sigma^2
    dinvgamma(sigmasq, shape = 1, rate = 1, log = TRUE)
  return(output)
}

s_target_fast <- function(theta, D, sigmasq, band.no) {
  output <- sbmds_ll_fast(theta, D, sigmasq, band.no) +
    # independent, standard Gaussian prior for theta
    sum(mvtnorm::dmvnorm(theta, log = TRUE)) +
    # inverse-gamma prior for sigma^2
    dinvgamma(sigmasq, shape = 1, rate = 1, log = TRUE)
  return(output)
}

##### delta function, rate of sd stepsize increase for adaptive MH/HMC

delta <- function(n) {
  return( min(0.01, n^(-0.5)) )
}

##### non-adaptive sigma, adaptive latent variable, Metropolis

sbmds_metropolis <- function(dims, maxIts, D, sigmasq, band.no,
                             targetAccept = 0.8, stepSize = 1) {

  # create the chain
  n <- dim(D)[1]
  chain <- array(0, dim = c(maxIts, n, dims))

  # specify the first random value
  chain[1, , ] <- mvtnorm::rmvnorm(n, mean = rep(0, dims), sigma = diag(dims))

  totalAccept <- rep(0, maxIts)
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 50   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)

  for (s in 2:maxIts) {

    thetaStar <- cbind(rnorm(n, mean = chain[s - 1, , 1], sd = stepSize),
                       rnorm(n, mean = chain[s - 1, , 2], sd = stepSize))
    # smaller proposal variance = more acceptances
    # proposal: normal dist
    u <- runif(1)
    # Metropolis, A = target(thetaStar)/target(previous iteration)
    # target on log scale
    logA <- s_target_no_sigma_prior(thetaStar, D, sigmasq, band.no) -
      s_target_no_sigma_prior(chain[s - 1, , ], D, sigmasq, band.no)

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

## fast

sbmds_metropolis_fast <- function(dims, maxIts, D, sigmasq, band.no,
                             targetAccept = 0.8, stepSize = 1) {

  # create the chain
  n <- dim(D)[1]
  chain <- array(0, dim = c(maxIts, n, dims))

  # specify the first random value
  chain[1, , ] <- mvtnorm::rmvnorm(n, mean = rep(0, dims), sigma = diag(dims))

  totalAccept <- rep(0, maxIts)
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 50   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)

  for (s in 2:maxIts) {

    thetaStar <- cbind(rnorm(n, mean = chain[s - 1, , 1], sd = stepSize),
                       rnorm(n, mean = chain[s - 1, , 2], sd = stepSize))
    # smaller proposal variance = more acceptances
    # proposal: normal dist
    u <- runif(1)
    # Metropolis, A = target(thetaStar)/target(previous iteration)
    # target on log scale
    logA <- s_target_no_sigma_prior_fast(thetaStar, D, sigmasq, band.no) -
      s_target_no_sigma_prior_fast(chain[s - 1, , ], D, sigmasq, band.no)

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

##### adaptive sigma & latent variable Metropolis

sbmds_mh <- function(dims, maxIts, D, band.no,
                     targetAccept = 0.8, targetAccept_Sigma = 0.8,
                     stepSize = 1, stepSizeSigma = 1) {

  # create the chain
  n <- dim(D)[1]
  chain <- array(0, dim = c(maxIts, n, dims)) # each iter is n x dims matrix
  sigma_chain <- rep(0, maxIts)

  # specify the first random value
  chain[1, , ] <- mvtnorm::rmvnorm(n, mean = rep(0, dims), sigma = diag(dims))
  sigma_chain[1] <- 1 # runif(1, min = 0.01, max = 0.99) # pick large number

  totalAccept <- rep(0, maxIts)
  totalAccept_Sigma <- rep(0, maxIts)
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  Acceptances_Sigma = 0
  SampBound = 50   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  SampCount_Sigma = 0

  for (s in 2:maxIts) {

    thetaStar <- cbind(rnorm(n, mean = chain[s - 1, , 1], sd = stepSize),
                       rnorm(n, mean = chain[s - 1, , 2], sd = stepSize))
    # smaller proposal variance = more acceptances
    # proposal: normal dist
    u <- runif(1)
    # Metropolis, A = target(thetaStar)/target(previous iteration)
    logA <- s_target(thetaStar, D, sigma_chain[s - 1], band.no) -
      s_target(chain[s - 1, , ], D, sigma_chain[s - 1], band.no)

    if(log(u) < logA) {
      chain[s, , ] <- thetaStar # ACCEPT !! # next iteration to thetastar
      totalAccept[s] <- 1
      Acceptances = Acceptances + 1 # acceptance counter
    } else {
      chain[s, , ] <- chain[s - 1, , ] # REJECT !!
      # next iteration stays at same iteration
    }

    # proposal distribution for sigma^2
    sigmaStar <- rtruncnorm(1, a = 0, b = Inf,
                            mean = sigma_chain[s - 1], sd = stepSizeSigma)

    # Metropolis, A2 = target(sigmaStar)/target(previous iteration)
    # comparing new and old sigma with current chain, no longer symmetric proposal
    logA2 <- s_target(chain[s, , ], D, sigmaStar, band.no) -
      s_target(chain[s, , ], D, sigma_chain[s - 1], band.no) +
      log(truncnorm::dtruncnorm(x = sigma_chain[s - 1], a = 0,
                                mean = sigmaStar, sd = stepSizeSigma)) -
      log(truncnorm::dtruncnorm(x = sigmaStar, a = 0,
                                mean = sigma_chain[s - 1], sd = stepSizeSigma))

    u2 <- runif(1)
    if(log(u2) < logA2) {
      sigma_chain[s] <- sigmaStar
      totalAccept_Sigma[s] <- 1
      Acceptances_Sigma <- Acceptances_Sigma + 1
    } else {
      sigma_chain[s] <- sigma_chain[s - 1]
    }

    SampCount <- SampCount + 1
    SampCount_Sigma <- SampCount_Sigma + 1

    # tune
    if (SampCount == SampBound) {
      AcceptRatio <- Acceptances / SampBound
      if ( AcceptRatio > targetAccept ) {
        stepSize <- stepSize * (1 + delta(s - 1)) # increase stepsize
      } else {
        stepSize <- stepSize * (1 - delta(s - 1)) # decrease stepsize
      }

      # reset Sampcount and Acceptances
      SampCount <- 0
      Acceptances <- 0
    }

    if (SampCount_Sigma == SampBound) {
      AcceptRatio_Sigma <- Acceptances_Sigma / SampBound
      if (AcceptRatio_Sigma > targetAccept_Sigma) {
        stepSizeSigma <- stepSizeSigma * (1 + delta(s - 1)) # increase stepsize
      } else {
        stepSizeSigma <- stepSizeSigma * (1 - delta(s - 1)) # decrease stepsize
      }

      # reset Sampcount and Acceptances
      SampCount_Sigma <- 0
      Acceptances_Sigma <- 0
    }

    if (s %% 100 == 0) cat("Iteration ", s, "\n","stepSize: ", stepSize, "\n",
                           "stepSizeSigma: ", stepSizeSigma, "\n")
  }

  cat("Acceptance rate: ", sum(totalAccept)/(maxIts - 1),
      "Acceptance rate for sigma: ", sum(totalAccept_Sigma)/(maxIts - 1))

  return(list(sigma_chain, chain))
}

## fast

sbmds_mh_fast <- function(dims, maxIts, D, band.no,
                     targetAccept = 0.8, targetAccept_Sigma = 0.8,
                     stepSize = 1, stepSizeSigma = 1) {

  # create the chain
  n <- dim(D)[1]
  chain <- array(0, dim = c(maxIts, n, dims)) # each iter is n x dims matrix
  sigma_chain <- rep(0, maxIts)

  # specify the first random value
  chain[1, , ] <- mvtnorm::rmvnorm(n, mean = rep(0, dims), sigma = diag(dims))
  sigma_chain[1] <- 1 # runif(1, min = 0.01, max = 0.99) # pick large number

  totalAccept <- rep(0, maxIts)
  totalAccept_Sigma <- rep(0, maxIts)
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  Acceptances_Sigma = 0
  SampBound = 50   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  SampCount_Sigma = 0

  for (s in 2:maxIts) {

    thetaStar <- cbind(rnorm(n, mean = chain[s - 1, , 1], sd = stepSize),
                       rnorm(n, mean = chain[s - 1, , 2], sd = stepSize))
    # smaller proposal variance = more acceptances
    # proposal: normal dist
    u <- runif(1)
    # Metropolis, A = target(thetaStar)/target(previous iteration)
    logA <- s_target_fast(thetaStar, D, sigma_chain[s - 1], band.no) -
      s_target_fast(chain[s - 1, , ], D, sigma_chain[s - 1], band.no)

    if(log(u) < logA) {
      chain[s, , ] <- thetaStar # ACCEPT !! # next iteration to thetastar
      totalAccept[s] <- 1
      Acceptances = Acceptances + 1 # acceptance counter
    } else {
      chain[s, , ] <- chain[s - 1, , ] # REJECT !!
      # next iteration stays at same iteration
    }

    # proposal distribution for sigma^2
    sigmaStar <- rtruncnorm(1, a = 0, b = Inf,
                            mean = sigma_chain[s - 1], sd = stepSizeSigma)

    # Metropolis, A2 = target(sigmaStar)/target(previous iteration)
    # comparing new and old sigma with current chain, no longer symmetric proposal
    logA2 <- s_target_fast(chain[s, , ], D, sigmaStar, band.no) -
      s_target_fast(chain[s, , ], D, sigma_chain[s - 1], band.no) +
      log(truncnorm::dtruncnorm(x = sigma_chain[s - 1], a = 0,
                                mean = sigmaStar, sd = stepSizeSigma)) -
      log(truncnorm::dtruncnorm(x = sigmaStar, a = 0,
                                mean = sigma_chain[s - 1], sd = stepSizeSigma))

    u2 <- runif(1)
    if(log(u2) < logA2) {
      sigma_chain[s] <- sigmaStar
      totalAccept_Sigma[s] <- 1
      Acceptances_Sigma <- Acceptances_Sigma + 1
    } else {
      sigma_chain[s] <- sigma_chain[s - 1]
    }

    SampCount <- SampCount + 1
    SampCount_Sigma <- SampCount_Sigma + 1

    # tune
    if (SampCount == SampBound) {
      AcceptRatio <- Acceptances / SampBound
      if ( AcceptRatio > targetAccept ) {
        stepSize <- stepSize * (1 + delta(s - 1)) # increase stepsize
      } else {
        stepSize <- stepSize * (1 - delta(s - 1)) # decrease stepsize
      }

      # reset Sampcount and Acceptances
      SampCount <- 0
      Acceptances <- 0
    }

    if (SampCount_Sigma == SampBound) {
      AcceptRatio_Sigma <- Acceptances_Sigma / SampBound
      if (AcceptRatio_Sigma > targetAccept_Sigma) {
        stepSizeSigma <- stepSizeSigma * (1 + delta(s - 1)) # increase stepsize
      } else {
        stepSizeSigma <- stepSizeSigma * (1 - delta(s - 1)) # decrease stepsize
      }

      # reset Sampcount and Acceptances
      SampCount_Sigma <- 0
      Acceptances_Sigma <- 0
    }

    if (s %% 100 == 0) cat("Iteration ", s, "\n","stepSize: ", stepSize, "\n",
                           "stepSizeSigma: ", stepSizeSigma, "\n")
  }

  cat("Acceptance rate: ", sum(totalAccept)/(maxIts - 1),
      "Acceptance rate for sigma: ", sum(totalAccept_Sigma)/(maxIts - 1))

  return(list(sigma_chain, chain))
}

##### non-adaptive sigma, adaptive latent variable HMC

sbmds_hmc_nonadapt_sigma <- function(dims, maxIts, D, sigmasq, band.no,
                                     targetAccept = 0.8, stepSize = 1) {

  n <- dim(D)[1]
  # initializations for latent variable
  chain <- array(0, dim = c(maxIts, n, dims))
  acceptances <- 0
  totalaccept <- rep(0, maxIts)
  SampCount <- 0
  SampBound <- 50   # current total samples before adapting radius

  L <- 20 # number of leapfrog steps

  # random starting point for latent variables and sigma
  chain[1, , ] <- mvtnorm::rmvnorm(n, mean = rep(0, dims), diag(dims))

  # U(q0) = - log posterior
  currentU <- - s_target_no_sigma_prior(chain[1, , ], D, sigmasq, band.no)

  for (i in 2:maxIts) {
    ####### update for latent variables - HMC
    proposalState <- chain[i - 1, , ] # q0
    momentum <- mvtnorm::rmvnorm(n, mean = rep(0, dims), diag(dims)) # p0
    currentK <- sum(momentum^2)/2 # valid bc independence; dimension K(p0)

    # leapfrog steps - obtain qt and pt
    momentum <- momentum + # half-step
      0.5 * stepSize * sbmds_grad(proposalState, D, sigmasq, band.no)

    for (l in 1:L) { # full step for p and q unless end of trajectory
      proposalState <- proposalState + stepSize * momentum # qt
      if (l != L) momentum <- momentum +
          stepSize * sbmds_grad(proposalState, D, sigmasq, band.no)
    }

    momentum <- momentum + # half-step
      0.5 * stepSize * sbmds_grad(proposalState, D, sigmasq, band.no)

    # quantities for accept/reject
    proposedU = - s_target_no_sigma_prior(proposalState, D, sigmasq, band.no) # U(qt)
    proposedK = sum(momentum^2)/2 # K(pt)
    u <- runif(1)

    if (log(u) < currentU - proposedU + currentK - proposedK) {
      chain[i, , ] <- proposalState # move pt to be p0 now
      currentU <- proposedU # update U(p0)
      totalaccept[i] <- 1
      acceptances <- acceptances + 1
    } else {
      chain[i, , ] <- chain[i - 1, , ] # keep p0
    }

    # adaptive stepsize for latent variables
    SampCount <- SampCount + 1

    if (SampCount == SampBound) {
      acceptRatio <- acceptances / SampBound
      if (acceptRatio > targetAccept) {
        stepSize <- stepSize * (1 + delta(i - 1)) # increase stepsize
      } else {
        stepSize <- stepSize * (1 - delta(i - 1)) # decrease stepsize
      }

      # reset Sampcount and Acceptances
      SampCount <- 0
      acceptances <- 0
    }

    if (i %% 100 == 0) cat("Iteration ", i, "\n","stepSize: ", stepSize)
  }

  cat("Acceptance rate: ", sum(totalaccept)/(maxIts - 1))

  return(chain)
}

##### adaptive sigma & latent variable HMC

sbmds_hmc <- function(dims, maxIts, D, band.no,
                      targetAccept = 0.8, targetAccept_Sigma = 0.8,
                      stepSize = 1, stepSizeSigma = 1) {

  n <- dim(D)[1]
  # initializations for latent variable
  chain <- array(0, dim = c(maxIts, n, dims))
  acceptances <- 0
  totalaccept <- rep(0, maxIts)
  SampCount <- 0
  SampBound <- 50   # current total samples before adapting radius

  # initializations for sigma
  sigma_chain <- rep(0, maxIts)
  totalaccept_Sigma <- rep(0, maxIts)
  acceptances_Sigma <- 0
  SampCount_Sigma <- 0

  L <- 20 # number of leapfrog steps

  # random starting point for latent variables and sigma
  chain[1, , ] <- mvtnorm::rmvnorm(n, mean = rep(0, dims), diag(dims))
  sigma_chain[1] <- 1
  # U(q0) = - log posterior
  currentU <- - s_target(chain[1, , ], D, sigma_chain[1], band.no)

  for (i in 2:maxIts) {
    ####### update for latent variables - HMC
    proposalState <- chain[i - 1, , ] # q0
    momentum <- mvtnorm::rmvnorm(n, mean = rep(0, dims), diag(dims)) # p0
    currentK <- sum(momentum^2)/2 # valid bc independence; dimension K(p0)

    # leapfrog steps - obtain qt and pt
    momentum <- momentum + # half-step
      0.5 * stepSize * sbmds_grad(proposalState, D, sigma_chain[i - 1], band.no)

    for (l in 1:L) { # full step for p and q unless end of trajectory
      proposalState <- proposalState + stepSize * momentum # qt
      if (l != L) momentum <- momentum +
          stepSize * sbmds_grad(proposalState, D, sigma_chain[i - 1], band.no)
    }

    momentum <- momentum + # half-step
      0.5 * stepSize * sbmds_grad(proposalState, D, sigma_chain[i - 1], band.no)

    # quantities for accept/reject
    proposedU = - s_target(proposalState, D, sigma_chain[i - 1], band.no) # U(qt)
    proposedK = sum(momentum^2)/2 # K(pt)
    u <- runif(1)

    if (log(u) < currentU - proposedU + currentK - proposedK) {
      chain[i, , ] <- proposalState # move pt to be p0 now
      currentU <- proposedU # update U(p0)
      totalaccept[i] <- 1
      acceptances <- acceptances + 1
    } else {
      chain[i, , ] <- chain[i - 1, , ] # keep p0
    }

    ####### update for sigma - adaptive Metropolis

    # proposal distribution for sigma
    sigmaStar <- rtruncnorm(1, a = 0, b = Inf, mean = sigma_chain[i - 1],
                            sd = stepSizeSigma)
    # comparing new and old sigma with current chain, not a symmetric proposal
    logA2 <- s_target(chain[i, , ], D, sigmaStar, band.no) + proposedU +
      # s_target(chain[i, , ], D, sigma_chain[i - 1], band.no) +
      log(truncnorm::dtruncnorm(x = sigma_chain[i - 1], a = 0,
                                mean = sigmaStar, sd = stepSizeSigma)) -
      log(truncnorm::dtruncnorm(x = sigmaStar, a = 0,
                                mean = sigma_chain[i - 1], sd = stepSizeSigma))

    u2 <- runif(1)
    if(log(u2) < logA2) {
      sigma_chain[i] <- sigmaStar
      totalaccept_Sigma[i] <- 1
      acceptances_Sigma <- acceptances_Sigma + 1
    } else {
      sigma_chain[i] <- sigma_chain[i - 1]
    }

    # adaptive stepsize for latent variables
    SampCount <- SampCount + 1

    if (SampCount == SampBound) {
      acceptRatio <- acceptances / SampBound
      if (acceptRatio > targetAccept) {
        stepSize <- stepSize * (1 + delta(i - 1)) # increase stepsize
      } else {
        stepSize <- stepSize * (1 - delta(i - 1)) # decrease stepsize
      }

      # reset Sampcount and Acceptances
      SampCount <- 0
      acceptances <- 0
    }

    # adaptive stepsize for sigma
    SampCount_Sigma <- SampCount_Sigma + 1

    if (SampCount_Sigma == SampBound) {
      acceptRatio_Sigma <- acceptances_Sigma / SampBound
      if (acceptRatio_Sigma > targetAccept_Sigma) {
        stepSizeSigma <- stepSizeSigma * (1 + delta(i - 1)) # increase stepsize
      } else {
        stepSizeSigma <- stepSizeSigma * (1 - delta(i - 1)) # decrease stepsize
      }

      # reset SampCount and Acceptances
      SampCount_Sigma <- 0
      acceptances_Sigma <- 0
    }

    if (i %% 100 == 0) cat("Iteration ", i, "\n","stepSize: ", stepSize, "\n",
                           "stepSizeSigma: ", stepSizeSigma, "\n")
  }

  cat("Acceptance rate: ", sum(totalaccept)/(maxIts - 1),
      "Acceptance rate sigma: ", sum(totalaccept_Sigma)/(maxIts - 1))

  return(list(sigma_chain, chain))
}

