### sBMDS example code

# import libraries
library(MassiveMDS)
library(mvtnorm)
library(truncnorm)
library(invgamma)

# set seeed
set.seed(12345)

# source functions
source("~/sparseBMDS/code/Rscripts/sBMDS_MassiveMDS_functions.R")

# inputs that can be varied
item <- 100 
band <- 10 # number of bands
sd_opt <- 0.1 # the amount of noise put on the data
dim <- 2

# simulate data
Z_ex <- rmvnorm(item, mean = rep(0, dim), sigma = diag(dim))
D_ex <- as.matrix(dist(Z_ex))

# add noise to distance matrix
noise <- matrix(rtruncnorm(item^2, mean = 0, sd = sd_opt), item, item) 
noise <- (noise + t(noise)) / 2
diag(noise) <- 0
D_noise <- noise + D_ex

# run simulation 
# with fixed sigma, i.e. set sigma to equal the amount of noise on the dist matrix
hmc_results1 <- sbmds_nonadapt_sigma_hmc(maxIts = 100000, dims = dim, 
                                         data = D_noise, precision = 1 / sd_opt^2, 
                                         bandwidth = band, targetAccept = 0.65, 
                                         stepSize = .01, thin = 100)

# without fixed sigma
hmc_results2 <- sbmds_hmc(maxIts = 100000, dims = dim, data = D_noise, 
                          bandwidth = band, targetAccept = 0.65, 
                          targetAccept_Sigma = 0.4, stepSize = .01, 
                          stepSizeSigma = 0.01, thin = 100)

# check for convergence
# NOTE: first value in thinned chain is from cmds estimate
# removed since want the first value at the thin"th" iteration

# traceplot of a latent variable
plot(hmc_results1[-1, 1, 1], type = "l")
plot(hmc_results2[[2]][-1, 1, 1], type = "l")

# traceplot for pairwise distances
hmcdist_mat1 <- array(0, dim = c(dim(hmc_results1)[1] - 1, item, item))
hmcdist_mat2 <- array(0, dim = c(dim(hmc_results2[[2]])[1] - 1, item, item))

for (j in 2:dim(hmc_results1)[1]){
  hmcdist_mat1[j - 1, , ] <- as.matrix(dist(hmc_results1[j, , ]))
  hmcdist_mat2[j - 1, , ] <- as.matrix(dist(hmc_results2[[2]][j, , ]))
}

par(mfrow=c(2,2))
for (i in 1:5){
  for (j in 10:15){ 
    plot(hmcdist_mat1[, i, j], type = "l")
    abline(h = D_ex[i, j], col = "red")
  }
}

par(mfrow=c(2,2))
for (i in 1:5){
  for (j in 10:15){ 
    plot(hmcdist_mat2[, i, j], type = "l")
    abline(h = D_ex[i, j], col = "green")
  }
}

# traceplot for sigma squared from adaptive HMC
plot(hmc_results2[[1]][-(1:1000)], type = "l") # chain not thinned
abline(h = sd_opt^2, col = "green")