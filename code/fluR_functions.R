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
#' SPEED-UPS: not writing to file as often, removing determinant calculations

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
dmatrixnorm <- function(X, Mu = NULL, Uinv, Vinv, gradient = FALSE) { #, logDetU, V) {
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
    #  (n / 2) * logDetV$modulus - (p / 2) * logDetU
    return(exponent) # remove attributes
  } # returns log-likelihood 
}

#### Potential function and gradient for HMC ####
# Takes MDS engine object, latent locations and other model parameters. Returns
# potential (proportional to log posterior) function or its gradient (if TRUE).
# log posterior = log-likelihood + log prior, 
# prior on locations: MN(Mu, treeVcv, traitVcv)
# prior on traitVcv: Wishart_p(traitT0^-1, d0)
# prior on precision: gamma(1, 1)

Potential <- function(engine, locations, treePrec, traitPrec, 
                      mdsprecision, gradient = FALSE) {
  
  engine <- MassiveMDS::updateLocations(engine, locations)
  engine <- MassiveMDS::setPrecision(engine, mdsprecision)
  
  if (gradient) {
    logPriorGrad <- dmatrixnorm(X = locations, Uinv = treePrec, 
                                #logDetU = treePrecDet, V = traitVcv, 
                                Vinv = traitPrec, gradient = gradient)
    logLikelihoodGrad <- MassiveMDS::getGradient(engine)
    return(-(logPriorGrad + logLikelihoodGrad)) # gradient of neg log-posterior
  }
  else {
    logPrior <- dmatrixnorm(X = locations, Uinv = treePrec, 
                            #logDetU = treePrecDet, V = traitVcv, 
                            Vinv = traitPrec) + # prior on X
      dgamma(mdsprecision, rate = 1, shape = 1, log = TRUE) # prior on precision
    logLikelihood <- MassiveMDS::getLogLikelihood(engine)
    return(-(logPrior + logLikelihood)) # negative log-posterior
  }
}

#### delta function for adaptive stepsize (HMC - latent variable)
delta <- function(n){return(min(0.01, n^(-0.5)))}

#### classical MDS for location start
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

#### MCMC (HMC) sampling to obtain latent locations ####
# returns list containing posterior samples, negative log likelihood values 
# and time to compute
hmcsampler <- function(n_iter,                  # number of MCMC samples
                       burnIn = 0,              # number of samples to throw
                       data,                    # distance matrix
                       beast = NULL,            # beast obj created by readbeast()
                       latentDimension = 2,     # dim of latent space
                       StepSize = 0.01,
                       StepSizeSigma = 0.01,
                       traitInvWeight = 1,      # how often to update TraitPrec
                       # how much tree to dictate covariance
                       priorRootSampleSize = 0.001, 
                       # if learnPrec = FALSE then mdsPrecision=1.47 (??)
                       mdsPrecision = 1,
                       sparse = FALSE,          # sparse gradient?
                       sparse_bands = 20,       # number of bands used in gradient          
                       truncation = TRUE,
                       targetAccept = 0.65,
                       targetAcceptSigma = 0.44,
                       threads = 1,
                       simd = 0, 
                       gpu = 0,
                       single = 0, 
                       thin = 100) {
  
  # Set up the parameters
  NumOfIterations = n_iter
  NumOfLeapfrog = 4
  P <- latentDimension
  N <- dim(beast$treeVcvRand)[1]
  #treePrec <- solve(beast$treeVcv)
  treeVcv.chol <- chol(beast$treeVcvRand)
  treePrec <- chol2inv(treeVcv.chol)
  #treePrecDet <- determinant(beast$treeVcv, logarithm = TRUE)$modulus
  beast$traitVcv <- diag(P) # initial est of traitVcv
  beast$d0 <- P # Wishart prior df on traitVcv
  beast$traitT0 <- diag(P) # Wishart prior cov on traitVcv
  
  # Allocate output space
  LocationSaved = list()
  locframe <- as.data.frame(matrix(NA, ncol = P * N)) # saving location
  Target = vector()
  mdsPrec <- rep(0, n_iter) # for the mds precision
  traitPrec <- array(0, dim = c(n_iter, P, P)) # for Vinv (p x p) from MN
  
  # Initial estimates
  set.seed(666)
  locations <- cmds(data, P) + matrix(rnorm(N * P, mean = 0, sd = .1), 
                                      nrow = N, ncol = P) #cmds(data, P)
  
  # Build reusable object to compute log-likelihood (gradient), i.e., engines
  # engine uses all bands, engine_lf uses sparse_bands
  engine <- engineStart(data, locations, N, P, bands = N, mdsPrecision, 
                        threads, simd, truncation, gpu, single)
  
  if (sparse == TRUE){
    engine_lf <- engineStart(data, locations, N, P, bands = sparse_bands, 
                             mdsPrecision, threads, simd, truncation, gpu, 
                             single)
  } else {
    engine_lf <- engine
  }
  
  mdsPrec[1] <- engine$precision
  traitPrec[1, , ] <- solve(beast$traitVcv) # beast$traitVcv = I_p
  
  # Control adaptive stepsizes
  Accepted = 0
  Proposed = 0
  AcceptPrec = 0
  ProposePrec = 0
  
  CurrentLocation = locations
  
  cat(paste0('Initial log-likelihood: ', MassiveMDS::getLogLikelihood(engine),
             '\n'))
  
  # Perform Hamiltonian Monte Carlo
  for (Iteration in 1:NumOfIterations) {
    
    # update CurrentU based on (new) location, precision, and traitPrec
    # U(q0) = -log posterior 
    CurrentU = Potential(engine = engine, locations = CurrentLocation, 
                         treePrec = treePrec,
                         #treePrecDet = treePrecDet, traitVcv = solve(traitPrec[1, , ]),
                         traitPrec = traitPrec[Iteration, , ], 
                         mdsprecision = mdsPrec[Iteration])
    
    if (Iteration > burnIn & Iteration %% thin == 0){
      
      traitCov <- solve(traitPrec[Iteration, , ])
      traitCovframe <- data.frame(r1c1 = traitCov[1, 1],
                                  r1c2 = traitCov[1, 2],
                                  r2c1 = traitCov[2, 1],
                                  r2c2 = traitCov[2, 2])
      write.table(traitCovframe, file = "latent_covL4D2.csv", append = TRUE,
                  row.names = FALSE, col.names = FALSE, sep = ',')
      
      chain_other <- data.frame(time = proc.time()[3], target = CurrentU,
                                mdsPrecision = mdsPrec[Iteration])
      write.table(chain_other, file = "chain_otherL4D2.csv", append = TRUE,
                  row.names = FALSE, col.names = FALSE, sep = ',')
      
      locframe[1, ] <- as.vector(CurrentLocation)
      
      write.table(as.data.frame(locframe), file = "latent_coordL4D2.csv",
                  append = TRUE, row.names = FALSE, col.names = FALSE,
                  sep = ',')
    }
    
    ProposedLocation = CurrentLocation
    
    # Sample the marginal momentum
    CurrentMomentum = matrix(rnorm(N * P), nrow = N, ncol = P) # p0
    ProposedMomentum = CurrentMomentum # K(p0)
    
    Proposed = Proposed + 1
    
    # half-step for momentum
    ProposedMomentum = ProposedMomentum -
      StepSize/2 * Potential(engine_lf, ProposedLocation, 
                             #treePrecDet = treePrecDet,
                             treePrec = treePrec,
                             #traitVcv = solve(traitPrec[Iteration, , ]),
                             traitPrec = traitPrec[Iteration, , ],
                             mdsprecision = mdsPrec[Iteration],
                             gradient = TRUE)
    
    # Simulate the Hamiltonian Dynamics
    for (StepNum in 1:NumOfLeapfrog) {
      # full step for position
      ProposedLocation = ProposedLocation + StepSize * ProposedMomentum
      if (StepNum != NumOfLeapfrog) {
        # full step for momentum until end of trajectory
        ProposedMomentum = ProposedMomentum -
          StepSize * Potential(engine_lf, ProposedLocation, 
                               #treePrecDet = treePrecDet,
                               treePrec = treePrec,
                               #traitVcv = solve(traitPrec[Iteration, , ]),
                               traitPrec = traitPrec[Iteration, , ],
                               mdsprecision = mdsPrec[Iteration],
                               gradient = TRUE)
      }
    }
    
    # half-step for momentum at the end
    ProposedMomentum = ProposedMomentum - 
      StepSize/2 * Potential(engine_lf, ProposedLocation, 
                             #treePrecDet = treePrecDet,
                             treePrec = treePrec,
                             #traitVcv = solve(traitPrec[Iteration, , ]),
                             traitPrec = traitPrec[Iteration, , ],
                             mdsprecision = mdsPrec[Iteration],
                             gradient = TRUE)
    
    ProposedMomentum = - ProposedMomentum # negate for symmetry
    
    # Compute the Potential
    ProposedU = Potential(engine, ProposedLocation, 
                          #treePrecDet = treePrecDet,
                          treePrec = treePrec,
                          #traitVcv = solve(traitPrec[Iteration, , ]),
                          traitPrec = traitPrec[Iteration, , ],
                          mdsprecision = mdsPrec[Iteration]) # U(qt)
    
    # Compute the Hamiltonian
    CurrentH = CurrentU + sum(CurrentMomentum^2) / 2 # = CurrentK K(q0)
    ProposedH = ProposedU + sum(ProposedMomentum^2) / 2 # = ProposedK K(qt)
    
    # Accept according to ratio
    Ratio = - ProposedH + CurrentH
    if (log(runif(1)) < Ratio){
      CurrentLocation = ProposedLocation
      CurrentU = ProposedU
      Accepted = Accepted + 1
    }
    
    # Save if sample is required
    if (Iteration > burnIn) {
      LocationSaved[[Iteration - burnIn]] = CurrentLocation
      Target[Iteration - burnIn] = CurrentU
    }
    
    # Show acceptance rate every 20 iterations
    if (Iteration %% 50 == 0) {
      cat(Iteration, "iterations completed. HMC acceptance rate: ",
          Accepted/Proposed,"\n", "StepSize", StepSize, "\n")
      
      if ((Accepted/Proposed) > targetAccept){
        StepSize <- StepSize * (1 + delta(Iteration))
      } else {
        StepSize <- StepSize * (1 - delta(Iteration))
      }
      
      Proposed = 0
      Accepted = 0
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
      prec_star <- truncnorm::rtruncnorm(1, a = 0, b = Inf,
                                         mean = mdsPrec[Iteration],
                                         sd = StepSizeSigma) # proposal dist, has to be pos
      ProposedU = Potential(engine, CurrentLocation, 
                            #treePrecDet = treePrecDet,
                            treePrec = treePrec,
                            #traitVcv = solve(traitPrec[Iteration, , ]),
                            traitPrec = traitPrec[Iteration, , ],
                            mdsprecision = prec_star) # g(prec_star)
      
      # comparing new and old sigma with current chain, not a symmetric proposal
      Ratio = - ProposedU + CurrentU +
        log(truncnorm::dtruncnorm(x = mdsPrec[Iteration], a = 0, b = Inf,
                                  mean = prec_star, sd = StepSizeSigma)) -
        log(truncnorm::dtruncnorm(x = prec_star, a = 0, b = Inf,
                                  mean = mdsPrec[Iteration], sd = StepSizeSigma))
      
      ProposePrec = ProposePrec + 1
      
      if (log(runif(1)) < Ratio) { 
        mdsPrec[Iteration + 1] <- prec_star
        AcceptPrec <- AcceptPrec + 1
      } else {
        mdsPrec[Iteration + 1] <- mdsPrec[Iteration]
        # engine <- MassiveMDS::setPrecision(engine, precision[Iteration])
        # replace engine back to previous precision, o.w. already updated above
      }
      
      if (Iteration %% 50 == 0) { # print MH acceptances
        cat(Iteration, "iterations completed. Prec acceptance rate: ",
            AcceptPrec/ProposePrec,"\n", "StepSizeSigma", StepSizeSigma, "\n")
        
        if ((AcceptPrec/ProposePrec) > targetAcceptSigma){
          StepSizeSigma <- StepSizeSigma * (1 + delta(Iteration))
        } else {
          StepSizeSigma <- StepSizeSigma * (1 - delta(Iteration))
        }
        
        ProposePrec = 0
        AcceptPrec = 0
      }
    }
    
    # Gibbs draw for traitPrec ~ Wishart(d0 + N, (traitT0 + t(X)%*%treePrec%*%X)^-1)
    if (Iteration < n_iter) {
      if (Iteration %% traitInvWeight == 0) { # same idea as thinning 
        traitPrec[Iteration + 1, , ] <- rWishart(1, beast$d0 + N, 
                                                 solve(beast$traitT0 + 
                                                         t(CurrentLocation) %*% 
                                                         treePrec %*% 
                                                         CurrentLocation))
      } else {
        traitPrec[Iteration + 1, , ] <- traitPrec[Iteration, , ]
      }
    }
  }
  
  time = proc.time() - timer
  acprat = dim(LocationSaved[!duplicated(LocationSaved)])[1] /
    (NumOfIterations - burnIn)
  
  return(list(samples = LocationSaved, target = Target, Time = time,
              acprat = acprat, precision = mdsPrec, traitPrec = traitPrec))
}