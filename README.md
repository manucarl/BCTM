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
All used data is freely available and made available through packages in the code.

## Code
- folder ./code/ contains helper functions used by BCTM in setting up the model
- folder code/nuts/ contains the NUTS code used for posterior inference. **Note that, currently, updates of the smoothing variances and results of the prior elicitation process for $\theta$ are hardcoded in many cases.**
- folder code/rcpp/ contains `rcpp` source files for posteriors (up to a constant) and gradients
- folder code/sims/ contains R code to reproduce all simulations and plots from the paper and the supplement
- folder code/framingham/ contains R code to reproduce all models and plots for the Framingham heart study illustration
- folder code/leukemia/ contains R code to reproduce all models and plots for the Leukemia lung cancer survival illustration
- folder code/veteran/ contains R code to reproduce all models and plots for the Veteran illustration
- folder code/mlt_sims/ contains R code from Hothorn et al. (2018) to simulate data

Used libraries are captured in the renv.lock file (after installing the `renv`package use `renv::init()`).

### Usage
All code in sims/framingham/leukemia/veteran is self-contained and can be used by simply running the scripts. Code can be used to reproduce the results in the paper and supplement, but is difficult to adapt for different tasks in its current state.


## Processed Files
Folder ./processed_data/ contains all results used in the paper except for files that are larger than GitHub's size limits.
