LandmarkLoglikelihood <- function(data, locations, precision, landmarks,
                                  truncation = TRUE, gradient = FALSE){
  locationCount <- dim(data)[1]
  sd <- 1 / sqrt(precision)
  logLikelihood <- 0
  gradlogLikelihood <- 0 * locations
  
  for (i in 1:landmarks) {
    for (j in i:locationCount) {
      if (i != j){
        mean <- as.numeric(dist(rbind(locations[i, ], locations[j, ])))
        if (gradient){
          if (truncation){
            # element in upper right triangle of landmark rows
            right_of_dii <- ((mean - data[i, j]) * precision + 
              dnorm(mean / sd) / (sd * pnorm(mean / sd))) * 
              (locations[i, ] - locations[j, ]) / mean
            # subtract when right of dii 
            gradlogLikelihood[i, ] <- colSums(rbind(gradlogLikelihood[i, ],
                                                    - right_of_dii))
            # add when left of dii (equal to locations[j, ] - locations[i, ])
            gradlogLikelihood[j, ] <- colSums(rbind(gradlogLikelihood[j, ], 
                                                    right_of_dii))
          } else {
            # element in upper right triangle of landmark rows
            right_of_dii <- ((mean - data[i, j]) * precision) * 
              (locations[i, ] - locations[j, ]) / mean
            # subtract when right of dii 
            gradlogLikelihood[i, ] <- colSums(rbind(gradlogLikelihood[i, ],
                                                    - right_of_dii))
            # add when left of dii (equal to locations[j, ] - locations[i, ])
            gradlogLikelihood[j, ] <- colSums(rbind(gradlogLikelihood[j, ], 
                                                    right_of_dii))
          }
        } else {
          # log-likelihood calc 
          logLikelihood <- sum(logLikelihood, dnorm(data[i, j], mean = mean, 
                                                    sd = sd, log = TRUE))
          if (truncation) {
            logLikelihood <- sum(logLikelihood, - log(pnorm(mean / sd)))
          }
        }
      }
    }
  }
  
  if (gradient) {
    return(gradlogLikelihood)
  } else {
    return(logLikelihood)
  }
}


set.seed(12345)
n <- 5
dims <- 2
Z.sim <- mvtnorm::rmvnorm(n, mean = rep(0, dims), sigma = diag(c(1, 1))) # dims = 2, n = 25
D.sim <- as.matrix(dist(Z.sim, diag = TRUE, upper = TRUE)) # n x n matrix of distances
noise <- matrix(rexp(n^2, rate = 1), n, n) # mean = 1 # increasing magnitude of noise
noise <- matrix(rtruncnorm(n^2, mean = 0, sd = .5), n, n)
noise <- (noise + t(noise)) / 2
diag(noise) <- 0
D.noise <- noise + D.sim

engine_ex <- MassiveMDS::createEngine(embeddingDimension = 2, locationCount = n,
                                      truncation = TRUE, tbb = 0, simd = 0, gpu = 0,
                                      single = 1, bandwidth = 1)

engine_ex <- MassiveMDS::setPairwiseData(engine_ex, D.sim)
engine_ex <- MassiveMDS::setPrecision(engine_ex, precision = 4)
engine_ex <- MassiveMDS::updateLocations(engine_ex, Z.sim)

MassiveMDS::getLogLikelihood(engine_ex)
MassiveMDS::getLogLikelihood(engine_ex, landmarks = T)
MassiveMDS::getGradient(engine_ex)
MassiveMDS::getGradient(engine_ex, landmarks = T)


MassiveMDS::computeLoglikelihood(D.sim, Z.sim, 4, 1, gradient = F) # R version loglike by bands
LandmarkLoglikelihood(data = D.sim, locations = Z.sim, precision = 4, 
                      landmarks = 3, gradient = T, truncation = T) # R version loglike by landmarks

landGrad <- function(latent) {
  LandmarkLoglikelihood(data = D.sim, locations = latent, precision = 4, 
                        landmarks = 3, truncation = T)
}

matrix(numDeriv::grad(landGrad, Z.sim), ncol = 2)
