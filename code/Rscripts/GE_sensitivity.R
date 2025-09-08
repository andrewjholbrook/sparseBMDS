# define arguments
args = commandArgs(trailingOnly=TRUE)
print(args)

items <- as.integer(args[1])
band <- as.double(args[2])
beta <- as.double(args[3])
iter <- as.double(args[4])
thinner <- as.double(args[5])
burnin <- as.double(args[6])

# import libraries
library(invgamma)
library(mvtnorm)
library(truncnorm)

# load necessary functions
source("code/Rscripts/GE_functions.R")

# set seed
set.seed(12345)

# simulate data
Z.sim <- mvtnorm::rmvnorm(items, mean = rep(0, 2), sigma = diag(2)) 
D.sim <- as.matrix(dist(Z.sim)) # n x n matrix of distances
# noise <- matrix(rtruncnorm(items^2, mean = 0, sd = sd_opt, a = 0), items, items) 
# noise <- (noise + t(noise)) / 2
# diag(noise) <- 0
b <- 0.2 / sqrt(2)
noise <- matrix(0, nrow = items, ncol = items)
# noise[upper.tri(noise)] <- rtruncnorm(items * (items - 1) / 2, mean = 0, 
#                                       sd = sd_opt, a = 0)
noise[upper.tri(noise)] <- VGAM::rlaplace(items * (items - 1) / 2, 
                                          location = 0, scale = b) # change scale 
noise <- (noise + t(noise))
D.noise <- noise + D.sim

# run simulation
sim <- bmds_hmc(maxIts = iter, dims = 2, data = D.noise, 
                bandwidth = band, beta = beta, landmark = TRUE,
                targetAccept = .65, targetAccept_Sigma = 0.44, 
                stepSize = .01, stepSizeSigma = .001, 
                thin = thinner, burnin = burnin)

# calculate mse and relative error for latent variables and sigma respectively
error <- c()

for (t in 1:(dim(sim$latent)[1])){
  if (items > 1000){
    subset <- sample(1:items, 1000)
    error[t] <- mean((as.matrix(dist(sim$latent[t, subset, ])) - 
                        as.matrix(dist(Z.sim[subset, ])))^2)
  } else {
    error[t] <- mean((as.matrix(dist(sim$latent[t, , ])) - D.sim)^2)
  }
}

error_results <- tibble(beta = beta, locations = items, band.no = band, 
                        method = "LM", mse = mean(error), 
                        rel_sigma = median(sim$sigma),
                        AR = sim$AR, ARS = sim$ARS)

write_csv(error_results, append = TRUE, file = "mse_GEmodel_V2.txt") #"~/sBMDS_results/mse_GEmodel.txt")
