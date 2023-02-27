# define arguments to be changed

args = commandArgs(trailingOnly=TRUE)
print(args)
sd_opt <- as.double(args[1])
band <- as.integer(args[2])
# oFile <- paste("errresults/", "sd_", sd_opt, "band#_", band, ".txt", sep="")

# import libraries
library(tidyverse)
library(invgamma)
library(mvtnorm)
library(truncnorm)
library(coda)
#library(gridExtra)

# load necessary functions
#source("sBMDS_functions.R")
source("R_scripts/sBMDS_functions.R")

### function for classic MDS

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

### function for Metropolis but using known latent variables as initial position

sbmds_metropolis_cheat <- function(dims, maxIts, D, sigmasq, band.no,
                                  targetAccept = 0.8, stepSize = 1) {

  # create the chain
  n <- dim(D)[1]
  chain <- array(0, dim = c(maxIts, n, dims))

  # specify the first random value
  #chain[1, , ] <- mvtnorm::rmvnorm(n, mean = rep(0, dims), sigma = diag(dims))
  #chain[1, , ] <- as.matrix(data_test2[ , -3])
  chain[1, , ] <- cmds(dist_test2, dims)

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

# set seed so that same dataset simulated for each sd
set.seed(12345)

# generate simulated data
l <- rep(c("blue", "green", "yellow", "orange", "red"), 5)

data_test <- tibble(x = rep(seq(-1, 1, by = .5), 5),
                    y = c(rep(-1, 5), rep(-.5, 5), rep(0, 5), rep(.5, 5), 
                          rep(1, 5)),
                    color = l)
dist_test <- as.matrix(dist(data_test[, 1:2], diag = TRUE, upper = TRUE))

## permuting items of the data
data_test2 <- data_test[sample(nrow(data_test)), ]
dist_test2 <- as.matrix(dist(data_test2[, 1:2], diag = TRUE, upper = TRUE))

## add noise to the data
size <- dim(dist_test2)[1]
noise <- matrix(rtruncnorm(size^2, mean = 0, sd = sd_opt), size, size)
noise <- (noise + t(noise)) / 2
diag(noise) <- 0
D.noise <- noise + dist_test2

error_grid <- function(sd, dims, size, maxIts, band.no, targetAccept, data,
                         stepSize, burnin){
  #data_sim <- array(0, dim = c(length(band.no), size, dims))
  #ees <- c()
  #noise <- matrix(rtruncnorm(size^2, mean = 0, sd = sd), size, size)
  #noise <- (noise + t(noise)) / 2
  #diag(noise) <- 0
  #D.noise <- noise + data
  # sim <- sbmds_metropolis_cheat(dims, maxIts, D = data,
  #                              sigmasq = sd^2, band.no, 
  #                              targetAccept, stepSize)
  
  sim <- sbmds_metropolis_cheat(dims, maxIts, D = data,
                   sigmasq = sd^2, band.no,
                   targetAccept, stepSize)
    
  ### calculate a distance matrix
  dist_mat_mh <- array(0, dim = c(maxIts, size, size))
    
  for (j in 1:maxIts){
    dist_mat_mh[j, , ] <- as.matrix(dist(sim[j, , ], diag = TRUE, upper = TRUE))
  }
  
  ess <- c()
  for (k in 1:size){
    ess[k] <- min(coda::effectiveSize(dist_mat_mh[-(1:burnin), -k, k]))
  }  
  #ees <- sum(coda::effectiveSize(dist_mat_mh[-(1:burnin), , 10]))
    ### look at minimum across all the (n 2) pairwise distances, > 100
  miness <- min(ess)
  #if (miness > 100){
  #  data_sim <- colMeans(sim[-(1:burnin), , ])
  #}
  data_sim <- colMeans(sim[-(1:burnin), , ])
  return(list(data_sim, miness))
  #return(list(ess, dist_mat_mh))
}

# results

error_results <- error_grid(sd = sd_opt, dims = 2, size = 25, 
                            maxIts = 300000, band.no = band, 
                            targetAccept = 0.238, data = D.noise, 
                            stepSize = .1, burnin = 1)

df_results <- tibble(sd = sd_opt, band.no = band, ees = error_results[[2]],
                     x.sim = error_results[[1]][ , 1],
                     y.sim = error_results[[1]][, 2],
                     color = data_test2$color)

#write.csv(df_results, file = oFile, sep = ",", row.names = FALSE)

write_csv(df_results, append = TRUE, file = "error_grid_results.txt")

# results checking
plot(error_results[[2]][-(1:50000), 25, 4], type = "l")
abline(h = dist_test2[25, 4], col = "red")
