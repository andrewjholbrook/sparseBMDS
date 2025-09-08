# define arguments
args = commandArgs(trailingOnly=TRUE)
print(args)

items <- as.integer(args[1])
HD_dim <- as.integer(args[2])
iter <- as.double(args[3])
thinner <- as.double(args[4])
burnin <- as.double(args[5])

# import libraries
library(MassiveMDS)
library(invgamma)
library(mvtnorm)
#library(tidyverse)
#library(dplyr)
library(truncnorm)
library(readr)
library(textmineR)

# load necessary functions
#source("~/sBMDS_scripts/sBMDS_MassiveMDS_functions.R")
source("Rscripts/sBMDS_MassiveMDS_functions.R")

# set seed
set.seed(12345)

# simulate data
Z.sim <- mvtnorm::rmvnorm(items, mean = rep(0, HD_dim), sigma = diag(HD_dim)) 
D.sim <- as.matrix(dist(Z.sim)) # n x n matrix of distances
noise <- matrix(0, nrow = items, ncol = items)
noise[upper.tri(noise)] <- rtruncnorm(items * (items - 1) / 2, 
                                      mean = 0, sd = .2, a = 0)
noise <- noise + t(noise)
D.noise <- noise + D.sim

# run full simulation
sim <- bmds_hmc(maxIts = iter, dims = HD_dim, data = D.noise, bandwidth = items,
                landmarks = FALSE, targetAccept = .65, 
                targetAccept_Sigma = 0.44, stepSize = .005, stepSizeSigma = .001, 
                thin = thinner, burnin = burnin)

# compute dist using latent locs
dist_array <- array(0, dim = c(dim(sim$latent)[1], items, items))
for (t in 1:(dim(sim$latent)[1])){
  dist_array[t, , ] <- as.matrix(dist(sim$latent[t, , ]))
}

# compute hellinger distances
band <- ceiling(HD_dim * sqrt(items))

# run simulation
sim2 <- bmds_hmc(maxIts = iter, dims = HD_dim, data = D.noise, bandwidth = band, 
                 landmarks = FALSE, targetAccept = .65, 
                 targetAccept_Sigma = 0.44, stepSize = .005, stepSizeSigma = .001, 
                 thin = thinner, burnin = burnin)
  
# compute dist using latent locs
dist_array2 <- array(0, dim = c(dim(sim2$latent)[1], items, items))
for (t in 1:(dim(sim2$latent)[1])){
    dist_array2[t, , ] <- as.matrix(dist(sim2$latent[t, , ]))
}
  
sumhellD <- 0
for (i in 1:(items - 1)){
  for (j in (i + 1):items){
    sumhellD <- sumhellD + CalcHellingerDist(dist_array[, i, j], dist_array2[, i, j])
  }
}
  
hellD <- sumhellD / (items * (items - 1) / 2)
print(hellD)
write_csv(data.frame(items = items, dims = HD_dim, band.no = band, 
                     hellD = hellD), append = TRUE, file = "BAND_HD.txt")
