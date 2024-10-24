# Sparse BMDS (sBMDS)

The folder `code` contains five different folders. We briefly describe each one and highlight important files:

1. `Rscript`:
   
`sBMDS_MassiveMDS_functions.R`: functions to compute sBMDS posterior distributions and MCMC algorithms using the `MassiveMDS` package. This R package efficiently evaluates the sBMDS log-likelihoods/gradients and integrates directly into BEAST. Because `MassiveMDS` contains C++ code, we found that it matches the theoretical speedups better than R code. We use this package in most of the Rscripts; however, users have the option to compute sBMDS log-likelihoods, gradients, posteriors, and MCMC algorithms directly in R in `sBMDS_functions.R`.

`MassiveMDS_rawllgrad.R`/ `llgrad_speedups.R`: output the computational time to evaluate the raw log-likelihood and gradient from the `MassiveMDS` package.

`MassiveMDS_mse.R`/ `LT_mse.R`/`HD_mse.R`: calculate mean MSE and median error for the latent variables and sigma respectively (under specified and misspecified cases)

`gaussian_grid.R`: simulate data from known Gaussian distributions and run HMC for a fixed sigma, returns the mean of the estimated latent locations across the iterations

`efficiency_comparison.R`: compare the time and the minimum effective sample size (ess) across different MCMC algorithms (MH, HMC) and BMDS variants (full BMDS, landmark sBMDS, banded sBMDS)

2. `data_graphs`: simulation output used for the plots in the sBMDS manuscript
   
3. `fluBMDS`: apply sBMDS variants to four influenza subtypes

`flufunc.R`: functions to compute the viral latent locations' priors and posteriors and run an HMC sampler

`fluSampler.R`: create the influenza subtype tree and run the HMC sampler

`fluPlots.R`: generate plots of viral diffusion rates and calculate Hellinger distances
   
4. `results_graphs`: the plots in the sBMDS manuscript
   
5. `rmd_files`: code used to generate the plots 



