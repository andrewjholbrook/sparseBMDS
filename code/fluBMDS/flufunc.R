#### Important parameters needed used within hmcsampler() ####
#' @param data N by N distance matrix.
#' @param locations N by P latent locations i.e. N P-dimensional latent locations.
#' @param N Number of locations and size of distance matrix.
#' @param P Dimension of latent locations.
#' @param precision MDS likelihood precision.
#' @param engine MDS engine object.
#' @param treeVcv N by N covariance.
#' @param traitVcv P by P covariance.
#' @param treePrec N by N precision.
#' @param traitPrec P by P precision.  
#' @param threads Number of CPU cores to be used.
#' @param simd For CPU implementation: no SIMD (\code{0}), SSE (\code{1}) or AVX (\code{2}).
#' @param truncation Likelihood includes truncation term? Defaults to \code{TRUE}.
#' @param gpu Which GPU to use? If only 1 available, use \code{gpu=1}. Defaults to \code{0}, no GPU.
#' @param single Set \code{single=1} if your GPU does not accommodate doubles.
#' 
#' 
#' seed up: not writing to file as often, speed up by removing determinants 

#### Initialize MDS engine object ####
engineStart <- function(data, locations, N, P, bands, precision = 1, 
                        threads, simd, truncation, gpu, single) {
  
  # Build reusable object
  engine <- MassiveMDS::createEngine(embeddingDimension = P, locationCount = N, 
                                     truncation = truncation, tbb = threads, 
                                     simd = simd, gpu = gpu, single = single,
                                     bandwidth = bands)
  # Set data once
  engine <- MassiveMDS::setPairwiseData(engine, as.matrix(data))
  # Call every time locations change
  engine <- MassiveMDS::updateLocations(engine, locations)
  # Call every time precision changes
  engine <- MassiveMDS::setPrecision(engine, precision = precision)
  return(engine)
}

#### Matrix Normal (MN) log-likelihood and gradient ####
# X ~ MN(Mu, U, V), U is n x n row covariance, V is p x p column covariance
dmatrixnorm <- function(X, Mu = NULL, Uinv, Vinv, gradient = FALSE) {
  #n <- dim(X)[1] # number of objects
  #p <- dim(X)[2] # latent dimension
  
  if (!is.null(Mu)) {
    X <- X - Mu # mean-centered
  } # else Mu = 0, and X - Mu = X
  
  if (gradient) {
    result <- - 0.5 * (t(Vinv %*% t(X) %*% Uinv) + Uinv %*% X %*% Vinv)
    return(unname(result)) # removes name attribute
  }
  else {
    product <- Vinv %*% t(X) %*% Uinv %*% X
    exponent <- - 0.5 * sum(diag(product)) # takes the trace
    #logDetV <- determinant(V, logarithm = TRUE)
    #logDetU <- determinant(U, logarithm = TRUE)
    #result <- exponent - (n * p / 2) * log(2 * pi) -
    #  (n / 2) * logDetV$modulus - (p / 2) * logDetU$modulus
    #return(as.vector(result)) # remove attributes
    return(exponent) # without normalizing constants
  } # returns log-likelihood 
}

#### Potential function and gradient for HMC ####
# Takes MDS engine object, latent locations and other model parameters. Returns
# potential (proportional to log posterior) function or its gradient (if TRUE).
# log posterior = log-likelihood + log prior, 
# prior on locations: MN(Mu, treeVcv, traitVcv)
# prior on traitPrec: Wishart_p(traitT0^-1, d0)
# prior on precision: gamma(1, 1)

Potential <- function(engine, locations, treePrec, traitPrec,
                      mdsprecision, gradient = FALSE, landmarks = FALSE) {
  
  engine <- MassiveMDS::updateLocations(engine, locations)
  engine <- MassiveMDS::setPrecision(engine, mdsprecision)
  
  if (gradient) {
    logPriorGrad <- dmatrixnorm(X = locations, Uinv = treePrec, 
                                Vinv = traitPrec, gradient = gradient)
    logLikelihoodGrad <- MassiveMDS::getGradient(engine, landmarks = landmarks)
    return(-(logPriorGrad + logLikelihoodGrad)) # gradient of neg log-posterior
  }
  else {
    logPrior <- dmatrixnorm(X = locations, Uinv = treePrec, 
                            Vinv = traitPrec) + # prior on X
      dgamma(mdsprecision, rate = 1, shape = 1, log = TRUE) # prior on precision
    logLikelihood <- MassiveMDS::getLogLikelihood(engine, landmarks = landmarks)
    return(-(logPrior + logLikelihood)) # negative log-posterior
  }
}

#### delta function for adaptive stepsize (HMC - latent variable)
delta <- function(n){return(min(0.01, n^(-0.5)))}

#### MCMC (HMC) sampling to obtain latent locations ####
# returns list containing posterior samples, negative log likelihood values 
# and time to compute
hmcsampler <- function(n_iter,                  # number of MCMC samples
                       burnIn = 0,              # number of samples to throw
                       data = data,             # distance matrix
                       virus = virus,           # type of virus
                       beast = beast,           # beast obj created by readbeast()
                       latentDimension = 2,     # dim of latent space
                       StepSize = 0.01,
                       StepSizeSigma = 0.01,
                       traitInvWeight = 1,      # how often to update TraitPrec
                       # how much tree to dictate covariance
                       priorRootSampleSize = 0.001,
                       mdsPrecision = 1,
                       sparse_bands = 20,       # number of bands for sparsity 
                       landmarks = FALSE,       # landmark or band method, default band
                       truncation = TRUE,
                       targetAccept = 0.65,
                       targetAcceptSigma = 0.45,
                       thin = 100,
                       threads = 0,
                       simd = 0, 
                       gpu = 0,
                       single = 0) {
  
  # Set up the parameters
  set.seed(666)
  NumOfLeapfrog = 20
  P <- latentDimension
  N <- dim(beast$treeVcvRand)[1]
  locations <- matrix(rnorm(N * P, 0, 1), nrow = N, ncol = P)
  treePrec <- solve(beast$treeVcvRand)   # Need inverse Vcovs
  #treeVcv.chol <- chol(beast$treeVcvRand)
  #treePrec <- chol2inv(treeVcv.chol)
  beast$traitVcv <- diag(P) * 5000 # initial est of traitVcv
  beast$d0 <- P # Wishart prior df on traitVcv
  beast$traitT0 <- diag(P) # Wishart prior cov on traitVcv
  
  # Allocate output space
  LocationSaved = list()
  locframe <- matrix(NA, ncol = P * N) # saving location
  Target = vector()
  mdsPrec <- rep(0, n_iter) # for the mds precision
  traitPrec <- array(0, dim = c(n_iter, P, P)) # for Vinv (p x p) from MN
  traitCovframe <- matrix(NA, ncol = P) # save dispersal rate
  thinCount <- 1 # thin target & locations
  
  # Build reusable object to compute log-likelihood (gradient)
  # initialize starting positions
  engine <- engineStart(data, locations, N, P, bands = sparse_bands, 
                        mdsPrecision, threads, simd, truncation, gpu, single)
  mdsPrec[1] <- engine$precision
  traitPrec[1, , ] <- solve(beast$traitVcv) # beast$traitVcv = I_p
  CurrentLocation = locations # q0
  
  cat(paste0('Initial log-likelihood: ', 
             MassiveMDS::getLogLikelihood(engine, landmarks = landmarks), '\n'))
  
  Accept = 0; Propose = 0;
  AcceptPrec = 0; ProposePrec = 0
  
  # Perform Hamiltonian Monte Carlo
  for (Iteration in 1:n_iter) {
    # update CurrentU based on (new) location, precision, and traitPrec, U(q0)
    CurrentU = Potential(engine = engine, locations = CurrentLocation, 
                         treePrec = treePrec, traitPrec = traitPrec[Iteration, , ], 
                         mdsprecision = mdsPrec[Iteration], landmarks = landmarks)
    # save current output to file
    if (Iteration > burnIn & Iteration %% thin == 0){

      # to be returned
      LocationSaved[[thinCount]] = CurrentLocation
      Target[thinCount] = CurrentU
      thinCount <- thinCount + 1
      
      traitCov <- solve(traitPrec[Iteration, , ]) 
      traitCovframe[1, ] <- diag(traitCov)
      write.table(as.data.frame(traitCovframe), 
                  file = paste0("latent_cov", virus, ".csv"),
                  append = TRUE, row.names = FALSE, 
                  col.names = FALSE, sep = ',')
      
      chain_other <- data.frame(time = proc.time()[3] - timer[3], 
                                target = CurrentU, 
                                mdsPrecision = mdsPrec[Iteration])
      write.table(chain_other, file = paste0("chain_other", virus, ".csv"),
                  append = TRUE,
                  row.names = FALSE, col.names = FALSE, sep = ',')
      
      locframe[1, ] <- as.vector(CurrentLocation)
      write.table(as.data.frame(locframe), 
                  file = paste0("latent_coord", virus, ".csv"),
                  append = TRUE, row.names = FALSE, col.names = FALSE,
                  sep = ',')
    }
    
    ProposedLocation = CurrentLocation
    
    # Sample the marginal momentum 
    CurrentMomentum = matrix(rnorm(N * P), nrow = N, ncol = P) # p0
    ProposedMomentum = CurrentMomentum # K(p0)
    
    Propose = Propose + 1
    
    # Simulate the Hamiltonian Dynamics
    for (StepNum in 1:NumOfLeapfrog) {
      # half-step for momentum
      ProposedMomentum = ProposedMomentum - 
        StepSize / 2 * Potential(engine = engine, locations = ProposedLocation, 
                                 treePrec = treePrec, 
                                 traitPrec = traitPrec[Iteration, , ], 
                                 mdsprecision = mdsPrec[Iteration], 
                                 gradient = TRUE, landmarks = landmarks)
      # full step for position
      ProposedLocation = ProposedLocation + StepSize * ProposedMomentum
      # half-step for momentum again
      ProposedMomentum = ProposedMomentum - 
        StepSize / 2 * Potential(engine = engine, locations = ProposedLocation, 
                                 treePrec = treePrec, 
                                 traitPrec = traitPrec[Iteration, , ], 
                                 mdsprecision = mdsPrec[Iteration],
                                 gradient = TRUE, landmarks = landmarks)
    }
    
    ProposedMomentum = - ProposedMomentum # negate for symmetry
    
    # Compute the Potential
    ProposedU = Potential(engine = engine, locations = ProposedLocation, 
                          treePrec = treePrec, 
                          traitPrec = traitPrec[Iteration, , ],
                          mdsprecision = mdsPrec[Iteration], 
                          landmarks = landmarks)
    
    # Compute the Hamiltonian
    CurrentH = CurrentU + sum(CurrentMomentum^2) / 2
    ProposedH = ProposedU + sum(ProposedMomentum^2) / 2
    
    # Accept according to ratio
    Ratio = - ProposedH + CurrentH
    if (log(runif(1)) < Ratio){
      CurrentLocation = ProposedLocation
      CurrentU = ProposedU
      Accept = Accept + 1
    }
    
    # Show acceptance rate every 20 iterations
    if (Iteration %% 50 == 0) {
      cat(Iteration, "iterations completed. HMC acceptance rate: ",
          Accept/Propose,"\n", "StepSize", StepSize, "\n")
      
      if ((Accept/Propose) > targetAccept){
        StepSize <- StepSize * (1 + delta(Iteration))
      } else {
        StepSize <- StepSize * (1 - delta(Iteration))
      }
      
      Propose = 0
      Accept = 0
    }
    
    # Start timer after burn-in
    if (Iteration == burnIn) { # If burnIn > 0
      cat("burn-in complete, now drawing samples ...\n")
      timer = proc.time()
    }
    if (burnIn == 0 & Iteration == 1) { # If burnIn = 0
      cat("burn-in complete, now drawing samples ...\n")
      timer = proc.time()
    }
    
    # MH step for residual precision
    if (Iteration < n_iter) {
      prec_star <- truncnorm:: rtruncnorm(1, a = 0, b = Inf, 
                                          mean = mdsPrec[Iteration], 
                                          sd = StepSizeSigma) # proposal dist, has to be pos
      ProposedU = Potential(engine = engine, locations = CurrentLocation, 
                            treePrec = treePrec, 
                            traitPrec = traitPrec[Iteration, , ], 
                            mdsprecision = prec_star, 
                            landmarks = landmarks) # g(prec_star)
      
      # comparing new and old sigma with current chain, not a symmetric proposal
      Ratio = - ProposedU + CurrentU + 
        log(truncnorm::dtruncnorm(x = mdsPrec[Iteration], a = 0, 
                                  mean = prec_star, sd = StepSizeSigma)) -
        log(truncnorm::dtruncnorm(x = prec_star, a = 0, 
                                  mean = mdsPrec[Iteration], sd = StepSizeSigma))
      
      # gamma prior for residual precision added to ratio
      ProposePrec <- ProposePrec + 1
      
      if (log(runif(1)) < Ratio) {
        mdsPrec[Iteration + 1] <- prec_star
        AcceptPrec <- AcceptPrec + 1
      } else {
        mdsPrec[Iteration + 1] <- mdsPrec[Iteration]
      }
      
      if (Iteration %% 50 == 0) { # print MH acceptances
        cat(Iteration, "iterations completed. Prec acceptance rate: ",
            AcceptPrec/ProposePrec,"\n", StepSizeSigma, "\n")
        
        if ((AcceptPrec/ProposePrec) > targetAcceptSigma){
          StepSizeSigma <- StepSizeSigma * (1 + delta(Iteration))
        } else {
          StepSizeSigma <- StepSizeSigma * (1 - delta(Iteration))
        }
        
        ProposePrec = 0
        AcceptPrec = 0
      }
    }
    
    # Gibbs draw for traitPrec if learnTraitPrec == TRUE
    # traitPrec ~ Wishart(d0 + N, (traitT0 + t(X)%*%treePrec%*%X)^-1)
    if (Iteration < n_iter) {
      if (Iteration %% traitInvWeight == 0) {
        traitPrec[Iteration + 1, , ] <- 
          rWishart(1, beast$d0 + N, 
                   solve(beast$traitT0 + 
                           t(CurrentLocation) %*% treePrec %*% CurrentLocation))
      } else {
        traitPrec[Iteration + 1, , ] <- traitPrec[Iteration,,]
      }
    }
    
  }
  
  time = proc.time() - timer
  acprat = dim(LocationSaved[!duplicated(LocationSaved)])[1] / 
    (n_iter - burnIn)
  
  return(list(samples = LocationSaved, target = Target, time = time, 
              acprat = acprat, precision = mdsPrec, 
              traitPrec = traitPrec))
}