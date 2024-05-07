# define arguments to be changed
args = commandArgs(trailingOnly=TRUE)
print(args)

items <- as.integer(args[1])
band <- as.double(args[2])
sd_opt <- as.double(args[3])

# import libraries
library(tidyverse)
library(invgamma)
library(mvtnorm)
library(truncnorm)
#library(gridExtra)

# load necessary functions
source("~/sBMDS_scripts/sBMDS_MassiveMDS_functions.R")

# simulate data
set.seed(12345)
latent <- rmvnorm(items, mean = rep(0, 2), sigma = diag(2))
D_latent <- as.matrix(dist(latent))

# add noise to data
noise <- matrix(0, nrow = items, ncol = items)
noise[upper.tri(noise)] <- rtruncnorm(items * (items - 1) / 2, mean = 0, 
                                      sd = sd_opt, a = 0)
noise <- noise + t(noise)
D.noise <- noise + D_latent

sim <- bmds_nonadapt_sigma_hmc(maxIts = 110000, dims = 2, data = D.noise, 
                               bandwidth = band, precision = (1 / sd_opt^2), 
                               landmarks = TRUE, targetAccept = 0.65, 
                               stepSize = 0.005, thin = 100, burnin = 10000)

# average across iterations
latent.mean <- colMeans(sim$latent[, , ])
    
# check ess aka proper convergence
if (items > 1000){
  size <- 1000
} else {
  size <- items
}

dist_sim <- array(0, dim = c(dim(sim$latent)[1], size, size))

for (t in 1:(dim(sim$latent)[1])){
  if (items > 1000){
    subset <- sample(1:item, 1000)
    dist_sim[t, , ] <- as.matrix(dist(sim$latent[t, subset, ]))
  } else {
    dist_sim[t, , ] <- as.matrix(dist(sim$latent[t, , ]))
  }
}

miness <- 110000
for (k in 1:(size - 1)){
  for (l in (k + 1):size){
    ess <- coda::effectiveSize(dist_sim[, k, l])
    if (ess < miness){
      miness <- ess
    }
  }
}

df_results <- tibble(sd = sd_opt, band.no = band, minees = miness,
                     x.sim = latent.mean[, 1],
                     y.sim = latent.mean[, 2])

#write.csv(df_results, file = oFile, sep = ",", row.names = FALSE)

write_csv(df_results, append = TRUE, file = "error_gaussian.txt")
