##### target functions
# (the log-posterior which is proportional to log-likelihood + log-priors)

s_target_no_sigma_prior <- function(locations, engine, landmarks) {
  # update engine to calculate likelihood at location of input
  engine <- MassiveMDS::updateLocations(engine, locations)
  # theta is the latent variable
  output <- MassiveMDS::getLogLikelihood(engine, landmarks) +
    # independent, Gaussian prior for theta centered at 0 & sd = 1
    sum(mvtnorm::dmvnorm(locations, mean = rep(0, engine$embeddingDimension),
                         sigma = diag(rep(1, engine$embeddingDimension)), 
                         log = TRUE))
  return(output)
}

s_target <- function(locations, engine, sigmasq, landmarks) {
  # update engine to calculate likelihood at location & sigmasq of input
  engine <- MassiveMDS::setPrecision(engine, precision = 1 / sigmasq)
  engine <- MassiveMDS::updateLocations(engine, locations)
  output <- MassiveMDS::getLogLikelihood(engine, landmarks) +
    # independent, Gaussian prior for theta centered at 0 & sd = 1
    sum(mvtnorm::dmvnorm(locations, mean = rep(0, engine$embeddingDimension),
                         sigma = diag(rep(1, engine$embeddingDimension)), 
                         log = TRUE)) +
    # inverse-gamma prior for sigma^2
    dinvgamma(sigmasq, shape = 1, rate = 1, log = TRUE)
  return(output)
}

##### gradient for target function

sbmds_grad <- function(locations, engine, landmarks){
  # update engine to calculate gradient at location of input
  engine <- MassiveMDS::updateLocations(engine, locations)
  gradient <- MassiveMDS::getGradient(engine, landmarks) - locations # from the MVN prior
  return(gradient)
}

sbmds_grad_sigma <- function(locations, engine, sigmasq, landmarks){
  # update engine to calculate gradient at location of input
  engine <- MassiveMDS::setPrecision(engine, precision = 1 / sigmasq)
  engine <- MassiveMDS::updateLocations(engine, locations)
  gradient <- MassiveMDS::getGradient(engine, landmarks) - locations # from the MVN prior
  return(gradient)
}

##### delta function, rate of sd stepsize for adaptive MH/HMC

delta <- function(n) {
  return(min(0.01, n^(-0.5)))
}

##### classical MDS function

cmds <- function(data, dims){
  # double center matrix
  dd <- data^2
  mndd <- sum(dd) / (nrow(dd)^2)
  rowdd <- dd*0 + rowMeans(dd)
  coldd <- t(dd*0 + colMeans(dd))
  B = - (dd - rowdd - coldd + mndd) / 2
  # decompose B by its eigenvalues and vectors
  eigendecomp <- eigen(B)
  # extract dims largest eigenvalues
  Lambda <- diag(eigendecomp$values[1:dims])
  # extract eigenvectors corresponding to eigenvalues
  E <- eigendecomp$vectors[ , 1:dims]
  # latent variable calculation, X = E %*% Lambda^(1/2)
  X <- E %*% sqrt(Lambda)
  return(X)
}

##### fixed sigma, adaptive Metropolis for latent variable

bmds_metropolis <- function(maxIts, dims, data, bandwidth, precision,
                            landmarks = FALSE, targetAccept = 0.8, 
                            stepSize = 1, thin = 1, burnin = 0) {
  n <- dim(data)[1]
  engine_test <- MassiveMDS::createEngine(embeddingDimension = dims, 
                                          locationCount = n, truncation = TRUE, 
                                          tbb = 0, simd = 0, gpu = 0,
                                          single = 0, bandwidth = bandwidth)
  engine_test <- MassiveMDS::setPairwiseData(engine_test, data)
  engine_test <- MassiveMDS::setPrecision(engine_test, precision)

  # create the chain
  latent_chain <- array(0, dim = c((maxIts - burnin) / thin, n, dims)) # n x dims matrices
  
  # specify the first random value
  # chain[1, , ] <- mvtnorm::rmvnorm(n, mean = rep(0, dims), sigma = diag(dims))
  loc <- cmds(data, dims) # classical MDS as initial value

  totalAccept <- rep(0, maxIts)
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 50   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  thinCount = 0 # to save thin samples

  for (s in 1:maxIts) {
    # proposal: normal dist
    thetaStar <- matrix(rnorm(n * dims, mean = loc, sd = stepSize),
                        ncol = dims, nrow = n)
    # smaller proposal variance = more acceptances
    # Metropolis, A = target(thetaStar)/target(previous iteration) 
    logA <- s_target_no_sigma_prior(locations = thetaStar, engine = engine_test, 
                                    landmarks = landmarks) -
      s_target_no_sigma_prior(locations = loc, engine = engine_test, 
                              landmarks = landmarks)
    
    u <- runif(1)
    
    if(log(u) < logA) {
      loc <- thetaStar # ACCEPT !! # next iteration to thetastar
      totalAccept[s] <- 1
      Acceptances = Acceptances + 1 # acceptance counter
    }  # else loc <- loc # REJECT !! i.e. next iteration stays at same iteration
    
    # adaptive stepsize
    SampCount <- SampCount + 1
    
    if (SampCount == SampBound) {
      AcceptRatio <- Acceptances / SampBound
      cat("Iteration ", s, "\n","stepSize: ", stepSize, "\n",
          "AR: ", AcceptRatio)
      if (AcceptRatio > targetAccept) {
        stepSize <- stepSize * (1 + delta(s - 1)) # increase stepsize
      } else {
        stepSize <- stepSize * (1 - delta(s - 1)) # decrease stepsize
      }

      # reset Sampcount and Acceptances
      SampCount <- 0
      Acceptances <- 0
    }
    
    if (s %% thin == 0 & s > burnin){
      thinCount <- thinCount + 1
      latent_chain[thinCount, , ] <- loc
    }
    
    if (s == burnin) { # If burnIn > 0
      cat("burn-in complete, now drawing samples ...\n")
      timer = proc.time()
    }
    if (burnin == 0 & s == 1) { # If burnIn = 0
      cat("burn-in complete, now drawing samples ...\n")
      timer = proc.time()
    }
  }
  
  time = proc.time() - timer
  cat("Acceptance rate: ", sum(totalAccept) / maxIts)
  return(list(latent = latent_chain, time = time))
}

##### fixed sigma, adaptive HMC with sparse gradient & likelihood for latent variables

bmds_nonadapt_sigma_hmc <- function(maxIts, dims, data, bandwidth, precision,
                                    landmarks = FALSE, targetAccept = 0.8, 
                                    stepSize = 1, thin = 1, burnin = 0) {
  n <- dim(data)[1]
  engine_test <- MassiveMDS::createEngine(embeddingDimension = dims, 
                                          locationCount = n, truncation = TRUE, 
                                          tbb = 0, simd = 0, gpu = 0,
                                          single = 0, bandwidth = bandwidth)
  engine_test <- MassiveMDS::setPairwiseData(engine_test, data)
  engine_test <- MassiveMDS::setPrecision(engine_test, precision)

  # initialization for latent variable
  latent_chain <- array(0, dim = c((maxIts - burnin) / thin, n, dims))
  target_chain <- c()
  acceptances <- 0
  totalaccept <- rep(0, maxIts)
  SampCount <- 0
  SampBound <- 50   # current total samples before adapting radius

  L <- 20 # number of leapfrog steps

  # random starting point for latent variables and sigma
  current_q <- cmds(data, dims) 
  thinCount <- 0 # allow for thinned chains

  # U(q0) = - log posterior
  currentU <- - s_target_no_sigma_prior(location = current_q, 
                                        engine = engine_test, 
                                        landmarks = landmarks)

  for (i in 1:maxIts) {
    
    ####### update for latent variables - HMC
    proposalState <- current_q # q0
    momentum <- mvtnorm::rmvnorm(n, mean = rep(0, dims), diag(dims)) # p0
    currentK <- sum(momentum^2)/2 # valid bc independence; K(p0)

    # leapfrog steps - obtain qt and pt
    momentum <- momentum + # half-step
      0.5 * stepSize * sbmds_grad(location = proposalState, 
                                  engine = engine_test, landmarks = landmarks)

    for (l in 1:L) { # full step for p and q unless end of trajectory
      proposalState <- proposalState + stepSize * momentum # qt
      if (l != L) momentum <- momentum +
          stepSize * sbmds_grad(location = proposalState, engine = engine_test,
                                landmarks = landmarks)
    }

    momentum <- momentum + # half-step
      0.5 * stepSize * sbmds_grad(location = proposalState, engine = engine_test,
                                  landmarks = landmarks)

    # quantities for accept/reject
    proposedU = - s_target_no_sigma_prior(location = proposalState, 
                                          engine = engine_test, 
                                          landmarks = landmarks) # U(qt)
    proposedK = sum(momentum^2)/2 # K(pt)
    
    u <- runif(1)
    if (log(u) < currentU - proposedU + currentK - proposedK) {
      current_q <- proposalState # move qt to be q0 now
      currentU <- proposedU # update U(p0)
      totalaccept[i] <- 1
      acceptances <- acceptances + 1
    } # else current_q and currentU stay the same

    # adaptive stepsize for latent variables
    SampCount <- SampCount + 1

    if (SampCount == SampBound) {
      acceptRatio <- acceptances / SampBound
      cat("Iteration ", i, "\n","stepSize: ", stepSize, "\n", 
          "AR: ", acceptRatio, "\n")
      if (acceptRatio > targetAccept) {
        stepSize <- stepSize * (1 + delta(i - 1)) # increase stepsize
      } else {
        stepSize <- stepSize * (1 - delta(i - 1)) # decrease stepsize
      }

      # reset Sampcount and Acceptances
      SampCount <- 0
      acceptances <- 0
    }
    
    if (i %% thin == 0 & i > burnin){
      thinCount <- thinCount + 1
      latent_chain[thinCount, , ] <- current_q
      target_chain[thinCount] <- currentU
    }
    
    # Start timer after burn-in
    if (i == burnin) { # If burnIn > 0
      cat("burn-in complete, now drawing samples ...\n")
      timer = proc.time()
    }
    if (burnin == 0 & i == 1) { # If burnIn = 0
      cat("burn-in complete, now drawing samples ...\n")
      timer = proc.time()
    }
  }
  
  time = proc.time() - timer

  cat("Acceptance rate: ", sum(totalaccept) / maxIts)
  return(list(latent = latent_chain, target = target_chain, time = time))
}

##### fixed sigma, adaptive HMC with sparse gradient & full likelihood for latent variables

bmds_nonadapt_sigma_sthmc <- function(maxIts, dims, data, bandwidth, precision, 
                                      targetAccept = 0.8, stepSize = 1, 
                                      thin = 1, burnin = 0) {
  n <- dim(data)[1]
  # create engine for sparse gradient
  engine_sparse <- MassiveMDS::createEngine(embeddingDimension = dims, 
                                            locationCount = n, 
                                            truncation = TRUE, tbb = 0, 
                                            simd = 0, gpu = 0, 
                                            single = 0, bandwidth = bandwidth)
  engine_sparse <- MassiveMDS::setPairwiseData(engine_sparse, data)
  engine_sparse <- MassiveMDS::setPrecision(engine_sparse, precision)
  
  # create engine for full likelihood
  engine_full <- MassiveMDS::createEngine(embeddingDimension = dims, 
                                          locationCount = n, truncation = TRUE, 
                                          tbb = 0, simd = 0, gpu = 0,
                                          single = 0, bandwidth = n)
  engine_full <- MassiveMDS::setPairwiseData(engine_full, data)
  engine_full <- MassiveMDS::setPrecision(engine_full, precision)

  # initialization for latent variable
  latent_chain <- array(0, dim = c((maxIts - burnin) / thin, n, dims))
  target_chain <- c()
  acceptances <- 0
  totalaccept <- rep(0, maxIts)
  SampCount <- 0
  SampBound <- 50   # current total samples before adapting radius

  L <- 20 # number of leapfrog steps

  # random starting point for latent variables and sigma
  current_q <- cmds(data, dims)
  thinCount <- 0 
  
  currentU <- - s_target_no_sigma_prior(location = current_q, 
                                        engine = engine_full)
  #timer = proc.time()
  for (i in 1:maxIts) {
    ####### update for latent variables - HMC
    proposalState <- current_q # q0
    momentum <- mvtnorm::rmvnorm(n, mean = rep(0, dims), diag(dims)) # p0
    currentK <- sum(momentum^2)/2 # valid bc independence; dimension K(p0)

    # leapfrog steps - obtain qt and pt
    momentum <- momentum + # half-step
      0.5 * stepSize * sbmds_grad(location = proposalState, engine = engine_sparse)

    for (l in 1:L) { # full step for p and q unless end of trajectory
      proposalState <- proposalState + stepSize * momentum # qt
      if (l != L) momentum <- momentum +
          stepSize * sbmds_grad(location = proposalState, engine = engine_sparse)
    }

    momentum <- momentum + # half-step
      0.5 * stepSize * sbmds_grad(location = proposalState, engine = engine_sparse)

    # quantities for accept/reject
    proposedU = - s_target_no_sigma_prior(location = proposalState, 
                                          engine = engine_full) # U(qt)
    proposedK = sum(momentum^2)/2 # K(pt)
    
    u <- runif(1)
    if (log(u) < currentU - proposedU + currentK - proposedK) {
      current_q <- proposalState # move pt to be p0 now
      currentU <- proposedU # update U(p0)
      totalaccept[i] <- 1
      acceptances <- acceptances + 1
    } # else current_q and currentU stay the same

    # adaptive stepsize for latent variables
    SampCount <- SampCount + 1

    if (SampCount == SampBound) {
      acceptRatio <- acceptances / SampBound
      cat("Iteration ", i, "\n","stepSize: ", stepSize, "\n", 
          "AR: ", acceptRatio, "\n")
      if (acceptRatio > targetAccept) {
        stepSize <- stepSize * (1 + delta(i - 1)) # increase stepsize
      } else {
        stepSize <- stepSize * (1 - delta(i - 1)) # decrease stepsize
      }

      # reset Sampcount and Acceptances
      SampCount <- 0
      acceptances <- 0
    }

    # save thinned values
    if (i %% thin == 0 & i > burnin){
      thinCount <- thinCount + 1
      latent_chain[thinCount, , ] <- current_q
      target_chain[thinCount] <- currentU
    }
  }

  #time = proc.time() - timer
  cat("Acceptance rate: ", sum(totalaccept) / maxIts)
  
  return(list(latent = latent_chain, target = target_chain, time = time))
}

##### adaptive sigma & latent variable Metropolis-Hastings

bmds_mh <- function(maxIts, dims, data, bandwidth, landmarks = FALSE, 
                    targetAccept = 0.8, targetAccept_Sigma = 0.8, 
                    stepSize = 1, stepSizeSigma = 1, thin = 1, burnin = 0) {
  
  n <- dim(data)[1]
  
  # initialization for latent variable
  latent_chain <- array(0, dim = c((maxIts - burnin) / thin, n, dims))
  acceptances <- 0
  totalaccept <- rep(0, maxIts)
  SampCount <- 0
  SampBound <- 50   # current total samples before adapting radius
  
  # initialization for sigma
  sigma_chain <- c()
  totalaccept_Sigma <- rep(0, maxIts)
  acceptances_Sigma <- 0
  SampCount_Sigma <- 0 # use the same sample bound
  thinCount <- 0
  
  # random starting point for latent variables and sigma
  loc <- cmds(data, dims) # first starting locations
  sigma2 <- .1 # first value of sigma^2
  
  # initialization of engine
  engine_test <- MassiveMDS::createEngine(embeddingDimension = dims, 
                                          locationCount = n, truncation = TRUE, 
                                          tbb = 0, simd = 0, gpu = 0,
                                          single = 0, bandwidth = bandwidth)
  engine_test <- MassiveMDS::setPairwiseData(engine_test, data)
  
  for (s in 1:maxIts) {
    # normal proposal for locations, smaller stepsize = more acceptances
    thetaStar <- matrix(rnorm(n * dims, mean = loc, sd = stepSize),
                        ncol = dims, nrow = n)
    # Metropolis, A = target(thetaStar)/target(previous iteration) 
    logA <- s_target(locations = thetaStar, engine = engine_test, 
                     sigmasq = sigma2, landmarks = landmarks) - 
      s_target(locations = loc, engine = engine_test, sigmasq = sigma2, 
               landmarks = landmarks)
    
    u <- runif(1)
    if(log(u) < logA) {
      loc <- thetaStar # next iteration to thetastar
      totalaccept[s] <- 1
      acceptances <- acceptances + 1 # acceptance counter
    } # else location stays the same 
    
    # proposal distribution for sigma^2
    sigmaStar <- rtruncnorm(1, a = 0, b = Inf,
                            mean = sigma2, sd = stepSizeSigma)

    # Metropolis, A2 = target(sigmaStar)/target(previous iteration)
    # comparing new and old sigma with current chain, no longer symmetric proposal
    logA2 <- s_target(locations = loc, engine = engine_test,
                      sigmasq = sigmaStar, landmarks = landmarks) -
      s_target(locations = loc, engine = engine_test, sigmasq = sigma2,
               landmarks = landmarks) +
      log(truncnorm::dtruncnorm(x = sigma2, a = 0, mean = sigmaStar,
                                sd = stepSizeSigma)) -
      log(truncnorm::dtruncnorm(x = sigmaStar, a = 0,
                                mean = sigma2, sd = stepSizeSigma))

    u2 <- runif(1)
    if(log(u2) < logA2) {
      sigma2 <- sigmaStar
      totalaccept_Sigma[s] <- 1
      acceptances_Sigma <- acceptances_Sigma + 1
    } # else sigma2 stays the same
    
    ####### update for sigma - approx Gibbs sampling
    # sigma^2 ~ IG(m/2 + a, SSR/2 + b) where prior sigma^2 ~ IG(a, b)
    # m <- n * (n - 1) / 2
    # m <- n * bandwidth - bandwidth * (bandwidth + 1) / 2
    #latent_dist <- as.matrix(dist(loc))
    # SSR <- sumr(loc, data, n)
    # SSR <- sum((latent_dist[upper.tri(latent_dist)] - data[upper.tri(data)])^2)
    # sigma2 <- invgamma::rinvgamma(1, shape = m / 2 + 1, rate = SSR / 2 + 1)
    
    SampCount <- SampCount + 1
    SampCount_Sigma <- SampCount_Sigma + 1
  
    if (SampCount == SampBound) {
      acceptRatio <- acceptances / SampBound
      cat("Iteration ", s, "\n","stepSize: ", stepSize, "\n",
          "AR: ", acceptRatio)
      if (acceptRatio > targetAccept) {
        stepSize <- stepSize * (1 + delta(s - 1)) # increase stepsize
      } else {
        stepSize <- stepSize * (1 - delta(s - 1)) # decrease stepsize
      }
      
      # reset Sampcount and Acceptances
      SampCount <- 0
      acceptances <- 0
    }
    
    if (SampCount_Sigma == SampBound) {
      acceptRatio_Sigma <- acceptances_Sigma / SampBound
      cat("stepSize: ", stepSizeSigma, "\n", "ARS: ", acceptRatio_Sigma)
      if (acceptRatio_Sigma > targetAccept_Sigma) {
        stepSizeSigma <- stepSizeSigma * (1 + delta(s - 1)) # increase stepsize
      } else {
        stepSizeSigma <- stepSizeSigma * (1 - delta(s - 1)) # decrease stepsize
      }

      # reset Sampcount and Acceptances
      SampCount_Sigma <- 0
      acceptances_Sigma <- 0
    }
    
    if (s %% thin == 0 & s > burnin){
      thinCount <- thinCount + 1
      latent_chain[thinCount, , ] <- loc
      sigma_chain[thinCount] <- sqrt(sigma2)
    }
  }
  
  cat("Acceptance rate: ", sum(totalaccept)/maxIts,
      "Acceptance rate sigma: ", sum(totalaccept_Sigma)/maxIts)
  return(list(latent = latent_chain, sigma = sigma_chain))
}

##### adaptive Metropolis for sigma
##### adaptive HMC with sparse gradient and likelihood for latent variables

bmds_hmc <- function(maxIts, dims, data, bandwidth, landmarks = FALSE,
                     targetAccept = 0.8, targetAccept_Sigma = 0.8, 
                     stepSize = 1, stepSizeSigma = 1, thin = 1, burnin = 0) {
  n <- dim(data)[1]
  # initialization of engine
  engine_test <- MassiveMDS::createEngine(embeddingDimension = dims, 
                                          locationCount = n, truncation = TRUE, 
                                          tbb = 0, simd = 0, gpu = 0,
                                          single = 0, bandwidth = bandwidth)
  engine_test <- MassiveMDS::setPairwiseData(engine_test, data)
  
  # engine_sig <- MassiveMDS::createEngine(embeddingDimension = dims,
  #                                        locationCount = n, truncation = TRUE,
  #                                        tbb = 0, simd = 0, gpu = 0,
  #                                        single = 0, bandwidth = bandwidth)
  # engine_sig <- MassiveMDS::setPairwiseData(engine_sig, data)
  
  # initialization for target 
  target_chain <- c()
  
  # initialization for latent variable
  latent_chain <- array(0, dim = c((maxIts - burnin) / thin, n, dims))
  acceptances <- 0
  totalaccept <- rep(0, maxIts)
  SampCount <- 0
  SampBound <- 50   # current total samples before adapting radius
  
  # initialization for sigma
  sigma_chain <- c()
  totalaccept_Sigma <- rep(0, maxIts)
  acceptances_Sigma <- 0
  SampCount_Sigma <- 0 # use the same sample bound
  
  # random starting point for latent variables and sigma
  current_q <- cmds(data, dims) # first starting locations
  current_s2 <- .1 # first value of sigma^2
  
  L <- 20 # number of leapfrog steps
  thinCount <- 0 # allow for thinning
  
  currentU <- - s_target(location = current_q, engine = engine_test, 
                         sigmasq = current_s2, landmarks = landmarks) 
  # U(q0) = - log posterior
  
  for (i in 1:maxIts) {
    ####### update for latent variables - HMC
    proposalState <- current_q # q0
    momentum <- mvtnorm::rmvnorm(n, mean = rep(0, dims), diag(dims)) # p0
    currentK <- sum(momentum^2)/2 # valid bc independence; dimension K(p0)
    
    # leapfrog steps - obtain qt and pt
    momentum <- momentum + # half-step
      0.5 * stepSize * sbmds_grad_sigma(location = proposalState, 
                                        engine = engine_test,
                                        sigmasq = current_s2, 
                                        landmarks = landmarks)
    
    for (l in 1:L) { # full step for p and q unless end of trajectory
      proposalState <- proposalState + stepSize * momentum # qt
      if (l != L) momentum <- momentum +
          stepSize * sbmds_grad_sigma(location = proposalState, 
                                      engine = engine_test,
                                      sigmasq = current_s2,
                                      landmarks = landmarks)
    }
    
    momentum <- momentum + # half-step
      0.5 * stepSize * sbmds_grad_sigma(location = proposalState, 
                                        engine = engine_test,
                                        sigmasq = current_s2,
                                        landmarks = landmarks)
    
    # quantities for accept/reject
    proposedU = - s_target(location = proposalState, engine = engine_test, 
                           sigmasq = current_s2, landmarks = landmarks) # U(qt)
    proposedK = sum(momentum^2)/2 # K(pt)
    
    u <- runif(1)
    if (log(u) < currentU - proposedU + currentK - proposedK) {
      current_q <- proposalState # move qt to be q0 now
      currentU <- proposedU # update U(p0)
      totalaccept[i] <- 1
      acceptances <- acceptances + 1
    } # else current_q and currentU remain the same
    
    ####### update for sigma - adaptive Metropolis
    
    # proposal distribution for sigma
    sigmaStar <- rtruncnorm(1, a = 0, b = Inf, mean = current_s2,
                            sd = stepSizeSigma)
    # comparing new and old sigma with current chain, not a symmetric proposal
    proposedU <- - s_target(location = current_q, engine = engine_test,
                            sigmasq = sigmaStar, landmarks = landmarks)

    logA2 <- -proposedU + currentU + 
      # s_target(location = current_q, engine = engine_sig,
      #                 sigmasq = sigmaStar, landmarks = landmarks) -
      # s_target(location = current_q, engine = engine_sig,
      #          sigmasq = current_s2, landmarks = landmarks) +
      log(truncnorm::dtruncnorm(x = current_s2, a = 0,
                                mean = sigmaStar, sd = stepSizeSigma)) -
      log(truncnorm::dtruncnorm(x = sigmaStar, a = 0,
                                mean = current_s2, sd = stepSizeSigma))

    u2 <- runif(1)
    if(log(u2) < logA2) {
      currentU <- proposedU
      # currentU <- - s_target(location = current_q, engine = engine_test,
      #                        sigmasq = sigmaStar, landmarks = landmarks)
      current_s2 <- sigmaStar
      totalaccept_Sigma[i] <- 1
      acceptances_Sigma <- acceptances_Sigma + 1
    } # current_s2 and currentU remain the same
    
    # adaptive stepsize for latent variables
    SampCount <- SampCount + 1
    
    if (SampCount == SampBound) {
      acceptRatio <- acceptances / SampBound
      cat("Iteration ", i, "\n","stepSize: ", stepSize, "\n",
          "AR: ", acceptRatio, "\n")
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
      cat("stepSizeSigma: ", stepSizeSigma, "\n",
          "ARS: ", acceptRatio_Sigma, "\n")
      if (acceptRatio_Sigma > targetAccept_Sigma) {
        stepSizeSigma <- stepSizeSigma * (1 + delta(i - 1)) # increase stepsize
      } else {
        stepSizeSigma <- stepSizeSigma * (1 - delta(i - 1)) # decrease stepsize
      }

      # reset SampCount and Acceptances
      SampCount_Sigma <- 0
      acceptances_Sigma <- 0
    }
    
    if (i %% thin == 0 & i > burnin){
      thinCount <- thinCount + 1
      latent_chain[thinCount, , ] <- current_q
      sigma_chain[thinCount] <- sqrt(current_s2)
      target_chain[thinCount] <- currentU
    }
  }
  
  cat("Acceptance rate: ", sum(totalaccept)/maxIts,
      "Acceptance rate sigma: ", sum(totalaccept_Sigma)/maxIts)
  
  return(list(latent = latent_chain, sigma = sigma_chain, target = target_chain,
              AR = sum(totalaccept) / maxIts,
              ARS = sum(totalaccept_Sigma) / maxIts))
}

##### adaptive HMC with sparse gradient and full likelihood for latent variables

bmds_sthmc <- function(maxIts, dims, data, bandwidth, landmarks = FALSE, 
                       targetAccept = 0.8, targetAccept_Sigma = 0.8, 
                       stepSize = 1, stepSizeSigma = 1, thin = 1, burnin = 0) {
  n <- dim(data)[1]
  
  # initialization for target 
  target_chain  <- c()
  
  # initialization for latent variable
  latent_chain <- array(0, dim = c((maxIts - burnin) / thin, n, dims))
  acceptances <- 0
  totalaccept <- rep(0, maxIts)
  SampCount <- 0
  SampBound <- 50   # current total samples before adapting radius
  
  # initialization for sigma
  sigma_chain <- c()
  totalaccept_Sigma <- rep(0, maxIts)
  acceptances_Sigma <- 0
  SampCount_Sigma <- 0 # use the same sample bound
  
  # random starting point for latent variables and sigma
  current_q <- cmds(data, dims)
  current_s2 <- .01 # first value of sigma^2
  
  # initialization of full engine for likelihood
  engine_full <- MassiveMDS::createEngine(embeddingDimension = dims, 
                                          locationCount = n, truncation = TRUE, 
                                          tbb = 0, simd = 0, gpu = 0,
                                          single = 0, bandwidth = n)
  engine_full <- MassiveMDS::setPairwiseData(engine_full, data)
  
  # initialization of sparse engine for gradient
  engine_sparse <- MassiveMDS::createEngine(embeddingDimension = dims, 
                                            locationCount = n, truncation = TRUE, 
                                            tbb = 0, simd = 0, gpu = 0,
                                            single = 0, bandwidth = bandwidth)
  engine_sparse <- MassiveMDS::setPairwiseData(engine_sparse, data)
  
  L <- 20 # number of leapfrog steps
  thin_count <- 0 # used to save thinned chains
  
  currentU <- - s_target(location = current_q, engine = engine_full, 
                         sigmasq = current_s2, landmarks = landmarks) 
  # U(q0) = - log posterior
  
  #timer = proc.time() # start timer at first iteration
  for (i in 1:maxIts) {
    
    ####### update for latent variables - HMC
    proposalState <- current_q # q0
    momentum <- mvtnorm::rmvnorm(n, mean = rep(0, dims), diag(dims)) # p0
    currentK <- sum(momentum^2)/2 # valid bc independence; dimension K(p0)
    
    # leapfrog steps - obtain qt and pt
    momentum <- momentum + # half-step
      0.5 * stepSize * sbmds_grad_sigma(location = proposalState, 
                                        engine = engine_sparse,
                                        sigmasq = current_s2,
                                        landmarks = landmarks)
    
    for (l in 1:L) { # full step for p and q unless end of trajectory
      proposalState <- proposalState + stepSize * momentum # qt, proposed q
      if (l != L) momentum <- momentum +
          stepSize * sbmds_grad_sigma(location = proposalState, 
                                      engine = engine_sparse,
                                      sigmasq = current_s2,
                                      landmarks = landmarks)
    }
    
    momentum <- momentum + # half-step
      0.5 * stepSize * sbmds_grad_sigma(location = proposalState, 
                                        engine = engine_sparse,
                                        sigmasq = current_s2,
                                        landmarks = landmarks)
    
    # negate trajectory to make proposal symmetric 
    momentum = -momentum # pt, proposed p
    
    # quantities for accept/reject
    proposedU = - s_target(location = proposalState, engine = engine_full, 
                           sigmasq = current_s2, landmarks = landmarks) # U(qt)
    proposedK = sum(momentum^2)/2 # K(pt)
    
    u <- runif(1)
    if (log(u) < currentU - proposedU + currentK - proposedK) {
      current_q <- proposalState # move qt to be q0 now
      currentU <- proposedU # update U(p0)
      totalaccept[i] <- 1
      acceptances <- acceptances + 1
    } # current_q and currentU stay the same
    
    ####### update for sigma - adaptive Metropolis-Hastings
    
    # proposal distribution for sigma
    sigmaStar <- rtruncnorm(1, a = 0, b = Inf, mean = current_s2,
                            sd = stepSizeSigma)
    # comparing new and old sigma with current chain, not a symmetric proposal
    proposedU <- - s_target(location = current_q, engine = engine_full, 
                            sigmasq = sigmaStar, landmarks = landmarks)
    
    logA2 <- - proposedU + currentU +
      log(truncnorm::dtruncnorm(x = current_s2, a = 0,
                                mean = sigmaStar, sd = stepSizeSigma)) -
      log(truncnorm::dtruncnorm(x = sigmaStar, a = 0,
                                mean = current_s2, sd = stepSizeSigma))
    
    u2 <- runif(1)
    if(log(u2) < logA2) {
      current_s2 <- sigmaStar
      currentU <- proposedU
      totalaccept_Sigma[i] <- 1
      acceptances_Sigma <- acceptances_Sigma + 1
    } # else current_s2 stays the same
    
    # adaptive stepsize for latent variables
    SampCount <- SampCount + 1
    
    if (SampCount == SampBound) {
      acceptRatio <- acceptances / SampBound
      cat("Iteration ", i, "\n","stepSize: ", stepSize, "\n", 
          "AR: ", acceptRatio, "\n")
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
      cat("stepSizeSigma: ", stepSizeSigma, "\n", 
          "ARS: ", acceptRatio_Sigma, "\n")
      if (acceptRatio_Sigma > targetAccept_Sigma) {
        stepSizeSigma <- stepSizeSigma * (1 + delta(i - 1)) # increase stepsize
      } else {
        stepSizeSigma <- stepSizeSigma * (1 - delta(i - 1)) # decrease stepsize
      }
      
      # reset SampCount and Acceptances
      SampCount_Sigma <- 0
      acceptances_Sigma <- 0
    }
    
    ## save output
    if (i %% thin == 0 & i > burnin){
      thin_count <- thin_count + 1
      latent_chain[thin_count, , ] <- current_q
      sigma_chain[thin_count] <- sqrt(current_s2)
      target_chain[thin_count] <- currentU
    }
  }
  
  #time = proc.time() - timer 
  cat("Acceptance rate: ", sum(totalaccept) / maxIts,
      "Acceptance rate sigma: ", sum(totalaccept_Sigma) / maxIts)
  
  return(list(latent = latent_chain, sigma = sigma_chain, target = target_chain))
}
