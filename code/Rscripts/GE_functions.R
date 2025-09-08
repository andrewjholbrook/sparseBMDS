### Likelihood function using bands
# beta = 1--> Gaussian dist; beta = 0 --> Laplace dist
## BANDED
banded_LL <- function(data, locations, sigmasq, bandwidth, 
                      beta, gradient = FALSE) {
  n <- nrow(locations)
  b <- sqrt(2 * sigmasq)
  if (bandwidth == n){
    stop("Error, bandwidth is max number of objects - 1.")
  }
  
  # precomputed once
  D_latent <- as.matrix(dist(locations, diag = TRUE, upper = TRUE))
  
  if (gradient) {
    grad <- matrix(0, nrow = n, ncol = ncol(locations))
    constant <- (1 + beta) / b^(1 + beta)
    
    for (j in 1:bandwidth) {
      i_idx <- 1:(n - j) # grabs entire band
      j_idx <- i_idx + j
      
      # vector differences
      x_diff <- locations[i_idx, , drop = F] - locations[j_idx, , drop = F]
      lat_dist <- D_latent[cbind(i_idx, j_idx)]
      delta_ij <- data[cbind(i_idx, j_idx)] - lat_dist
      weight <- constant * (abs(delta_ij)^beta) * sign(delta_ij) / lat_dist
      grad[i_idx, ] <- grad[i_idx, ] + weight * x_diff
      grad[j_idx, ] <- grad[j_idx, ] - weight * x_diff
    }
    return(grad)
    
  } else {
    ssr_sum <- 0
    for (j in 1:bandwidth) {
      i_idx <- 1:(n - j) # grabs entire band
      j_idx <- i_idx + j 
      lat_dist <- D_latent[cbind(i_idx, j_idx)]
      obs_dist <- data[cbind(i_idx, j_idx)]
      ssr <- (abs(obs_dist - lat_dist) / b)^(1 + beta)
      ssr_sum <- ssr_sum + sum(ssr, na.rm = TRUE)
    }
    m <- n * bandwidth - bandwidth * (bandwidth + 1) / 2
    return(-m * log(b) - ssr_sum)
  }
}

## LANDMARK
landmark_LL <- function(data, locations, sigmasq, bandwidth, 
                        beta, gradient = FALSE) {
  n <- nrow(locations)
  b <- sqrt(2 * sigmasq)
  
  # precomputed once
  D_latent <- as.matrix(dist(locations, diag = TRUE, upper = TRUE))
  
  if (gradient) {
    grad <- matrix(0, nrow = n, ncol = ncol(locations))
    constant <- (1 + beta) / b^(1 + beta)
    
    for (j in 1:bandwidth) {
      i_idx <- rep(j, n - j) # grabs row
      j_idx <- (j + 1):n
      
      # vector differences
      x_diff <- locations[i_idx, , drop = F] - locations[j_idx, , drop = F]
      lat_dist <- D_latent[cbind(i_idx, j_idx)]
      delta_ij <- data[cbind(i_idx, j_idx)] - lat_dist
      weight <- constant * (abs(delta_ij)^beta) * sign(delta_ij) / lat_dist
      grad[j, ] <- grad[j, ] + colSums(weight * x_diff)
      grad[j_idx, ] <- grad[j_idx, ] - weight * x_diff
    }
    return(grad)
    
  } else {
    ssr_sum <- 0
    for (j in 1:bandwidth) {
      i_idx <- rep(j, n - j) # grabs row
      j_idx <- (j + 1):n
      lat_dist <- D_latent[cbind(i_idx, j_idx)]
      obs_dist <- data[cbind(i_idx, j_idx)]
      ssr <- (abs(obs_dist - lat_dist) / b)^(1 + beta)
      ssr_sum <- ssr_sum + sum(ssr, na.rm = TRUE)
    }
    m <- n * bandwidth - bandwidth * (bandwidth + 1) / 2
    return(-m * log(b) - ssr_sum)
  }
}

#### target function
s_target <- function(data, locations, sigmasq, bandwidth, beta, landmark) {
  # calculate likelihood 
  dims <- dim(locations)[2]
  if (landmark) {
    loglik <- landmark_LL(data, locations, sigmasq, bandwidth, beta)
  } else {
    loglik <- banded_LL(data, locations, sigmasq, bandwidth, beta)
  }
  output <- loglik +
    # independent, Gaussian prior for theta centered at 0 & sd = 1
    sum(mvtnorm::dmvnorm(locations, mean = rep(0, dims),
                         sigma = diag(rep(1, dims)),
                         log = TRUE)) +
    # inverse-gamma prior for sigma^2
    invgamma::dinvgamma(sigmasq, shape = 1, rate = 1, log = TRUE)
  return(output)
}

#### gradient for target function
sbmds_grad_sigma <- function(data, locations, sigmasq, bandwidth, beta, 
                             landmark) {
  # calculate gradient at location of input
  if (landmark) {
    gradient <- landmark_LL(data, locations, sigmasq, bandwidth, beta,
                            gradient = TRUE) - locations # from MVN prior
  } else {
    gradient <- banded_LL(data, locations, sigmasq, bandwidth, beta, 
                          gradient = TRUE) - locations # from MVN prior
  }
  return(gradient)
}

##### delta function, rate of sd stepsize increase for adaptive MH/HMC
delta <- function(n) {
  return( min(0.01, n^(-0.5)) )
}

##### classical MDS function
cmds <- function(D, dims){
  # double center matrix
  dd <- D^2
  mndd <- mean(dd)
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

#### adaptive HMC (sigma & latent locations)
bmds_hmc <- function(maxIts, dims, data, bandwidth, beta, landmark = FALSE,
                     targetAccept = 0.8, targetAccept_Sigma = 0.8, 
                     stepSize = 1, stepSizeSigma = 1, thin = 1, burnin = 0) {
  n <- dim(data)[1]
  
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
  SampCount_Sigma <- 0
  
  # random starting point for latent variables and sigma
  current_q <- cmds(data, dims) + 
    matrix(rnorm(n * dims, mean = 0, sd = 0.01), nrow = n, ncol = dims)
  # first starting locations
  current_s2 <- .1 # first value of sigma^2
  # U(q0) = - log posterior
  currentU <- - s_target(data = data, locations = current_q, 
                         sigmasq = current_s2, bandwidth = bandwidth, 
                         beta = beta, landmark = landmark) 
  
  L <- 20 # number of leapfrog steps
  thinCount <- 0 # allow for thinning
  
  for (i in 1:maxIts) {
    ####### update for latent variables - HMC
    proposalState <- current_q # q0
    momentum <- mvtnorm::rmvnorm(n, mean = rep(0, dims), diag(dims)) # p0
    currentK <- sum(momentum^2)/2 # valid bc independence; dimension K(p0)
    
    # leapfrog steps - obtain qt and pt
    momentum <- momentum + # half-step
      0.5 * stepSize * sbmds_grad_sigma(data = data, locations = proposalState, 
                                        sigmasq = current_s2, 
                                        bandwidth = bandwidth, beta = beta, 
                                        landmark = landmark)
    
    for (l in 1:L) { # full step for p and q unless end of trajectory
      proposalState <- proposalState + stepSize * momentum # qt
      if (l != L) momentum <- momentum +
          stepSize * sbmds_grad_sigma(data = data, locations = proposalState, 
                                      sigmasq = current_s2, 
                                      bandwidth = bandwidth, beta = beta, 
                                      landmark = landmark)
    }
    
    momentum <- momentum + # half-step
      0.5 * stepSize * sbmds_grad_sigma(data = data, locations = proposalState, 
                                        sigmasq = current_s2, 
                                        bandwidth = bandwidth, beta = beta, 
                                        landmark = landmark)
    
    # quantities for accept/reject
    proposedU = - s_target(data = data, locations = proposalState, 
                           sigmasq = current_s2, bandwidth = bandwidth, 
                           beta = beta, landmark = landmark) # U(qt)
    proposedK = sum(momentum^2)/2 # K(pt)
    
    u <- runif(1)
    if (log(u) < currentU - proposedU + currentK - proposedK) {
      current_q <- proposalState # move pt to be p0 now
      currentU <- proposedU # update U(p0)
      totalaccept[i] <- 1
      acceptances <- acceptances + 1
    } # else current_q and currentU remain the same
    
    ####### update for sigma - adaptive Metropolis
    
    # proposal distribution for sigma
    sigmaStar <- rtruncnorm(1, a = 0, b = Inf, mean = current_s2,
                            sd = stepSizeSigma)
    # comparing new and old sigma with current chain, not a symmetric proposal
    proposedU <- - s_target(data = data, location = current_q,
                            sigmasq = sigmaStar, 
                            bandwidth = bandwidth, beta = beta, 
                            landmark = landmark)
    
    logA2 <- - proposedU + currentU + 
      log(truncnorm::dtruncnorm(x = current_s2, a = 0,
                                mean = sigmaStar, sd = stepSizeSigma)) -
      log(truncnorm::dtruncnorm(x = sigmaStar, a = 0,
                                mean = current_s2, sd = stepSizeSigma))
    
    u2 <- runif(1)
    if(log(u2) < logA2) {
      currentU <- proposedU
      current_s2 <- sigmaStar
      totalaccept_Sigma[i] <- 1
      acceptances_Sigma <- acceptances_Sigma + 1
    } # current_s2 and currentU stay the same
    
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
  
  cat("Acceptance rate: ", sum(totalaccept)/(maxIts - 1),
      "Acceptance rate sigma: ", sum(totalaccept_Sigma)/(maxIts - 1))
  
  return(list(latent = latent_chain, sigma = sigma_chain, target = target_chain,
              AR = sum(totalaccept) / maxIts,
              ARS = sum(totalaccept_Sigma) / maxIts))
}