Bayesian Conditional Transformation Models
================

## Abstract
Recent developments in statistical regression methodology shift away from pure mean
regression towards distributional regression models. One important strand thereof is
that of conditional transformation models (CTMs). CTMs infer the entire conditional
distribution directly by applying a transformation function to the response condi-
tionally on a set of covariates towards a simple log-concave reference distribution.
Thereby, CTMs allow not only variance, kurtosis or skewness but the complete con-
ditional distribution to depend on the explanatory variables. We propose a Bayesian
notion of conditional transformation models (BCTMs) focusing on exactly observed
continuous responses, but also incorporating extensions to randomly censored and
discrete responses. Rather than relying on Bernstein polynomials that have been
considered in likelihood-based CTMs, we implement a spline-based parametrization
for monotonic effects that are supplemented with smoothness priors. Furthermore,
we are able to benefit from the Bayesian paradigm via easily obtainable credible in-
tervals and other quantities without relying on large sample approximations. A simu-
lation study demonstrates the competitiveness of our approach against its likelihood-
based counterpart but also Bayesian additive models of location, scale and shape and
Bayesian quantile regression. Two applications illustrate the versatility of BCTMs
in problems involving real world data, again including the comparison with various
types of competitors.

## Data
All used data is freely available through packages:

- Framingham illustration: data("Cholesterol", package="qrLMM")
- Leukemia survival illustration: data("LeukSurv", package="spBayesSurv")
- Veteran lung cancer illustration: data(cancer, package="survival")

and via csv in the folder ./data.


## Code
- folder ./code/ contains helper functions used by BCTM in setting up the model
- folder code/nuts/ contains the NUTS code used for posterior inference. **Note that, currently, updates of the smoothing variances and results of the prior elicitation process for $\theta$ are hardcoded in many cases.**
- folder code/rcpp/ contains `rcpp` source files for posteriors (up to a constant) and gradients
- folder code/sims/ contains R code to reproduce all simulations and plots from the paper and the supplement
- folder code/framingham/ contains R code to reproduce all models and plots for the Framingham heart study illustration
- folder code/leukemia/ contains R code to reproduce all models and plots for the Leukemia survival illustration
- folder code/veteran/ contains R code to reproduce all models and plots for the Veteran illustration
- folder code/mlt_sims/ contains R code from Hothorn et al. (2018) to simulate data

All scripts include a helper function that installs and loads used libraries. A majority of (but not all) libraries are captured in the renv.lock file (after installing the `renv`package use `renv::init()`).

### Usage
All code in sims/framingham/leukemia/veteran is self-contained and can be used by simply running the scripts. Code can be used to reproduce the results in the paper and supplement, but is difficult to adapt for different tasks in its current state. The different files in code/nuts/ differ only in the updates of the smoothing variances and prior elicitation schemes (nuts.R is the default used in most of the applications).
Naming:

- _sd refers to scale-dependent prior
- _te refers to tensor spline
- _re refers to random effect
- _vcm refers to varying coefficient model

### Dependencies
- splines (4.21)
- dplyr (1.0.10)
- tidyverse (1.3.2)
- mlt (1.4-2)
- tram (0.7-2)
- tramME (1.0.3)
- scam (1.2-13)
- mgcv (1.8-41)
- bamlss (1.1-8)
- Rcpp (1.0.9) RcppArmadillo (0.11.2.3.1) RcppEigen (0.3.3.92)
- qrLMM (2.1)
- viridis (0.6.2)
- MASS (7.3-58.1)
- mvtnorm (1.1-3)
- cowplot  (1.1.1)
- gamboostLSS (2.0-6)
- lattice (0.20-45)
- microbenchmark (1.4.9)
- progress (1.2.2)
- devtools (2.4.4)
- doParallel (1.0.17)
- survival (3.4-0)
- spBayesSurv (1.1.6)
- Matrix (1.5-1)
- MCMCpack (1.6-3)
- profvis (0.3.7)
- sf (1.0-8)
- rgeos (0.5-9)
- RhpcBLASctl (0.21-247.1)
- Loo (2.5.1)
- BayesX (0.3-1.1)
- caret (6.0-93)
- sdPrior (1.0-0)   
- scoringutils (1.0.1)


## Processed Files
Folder ./processed_data/ contains all results used in the paper except for files that are larger than GitHub's size limits (due to inclusion of all samples).
