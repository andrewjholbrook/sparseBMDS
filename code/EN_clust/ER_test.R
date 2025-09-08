# define arguments
args = commandArgs(trailingOnly=TRUE)
print(args)

band <- as.double(args[1])
iter <- as.double(args[2])
thinner <- as.double(args[3])
burnin <- as.double(args[4])

# import libraries
library(MassiveMDS)
library(invgamma)
library(mvtnorm)
#library(tidyverse)
#library(dplyr)
library(truncnorm)
library(readr)
#library(textmineR)

# load necessary functions
#source("~/sBMDS_scripts/sBMDS_MassiveMDS_functions.R")
source("Rscripts/sBMDS_MassiveMDS_functions.R")

# read in ArXiv data
arXiv_data <- as.matrix(read_csv("dist_mini.csv", col_names = FALSE))
N <- nrow(arXiv_data)
#n <- 50

# set seed
set.seed(12345)

# run simultation
sim <- bmds_hmc(maxIts = iter, dims = 2, data = arXiv_data, bandwidth = band,
                landmarks = FALSE, targetAccept = .65, 
                targetAccept_Sigma = 0.44, stepSize = .01, stepSizeSigma = .01, 
                thin = thinner, burnin = burnin)