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
- folder ./nuts/ contains the NUTS code used for posterior inference. **Note that, currently, updates of the smoothing variances and results of the prior elicitation process for $\theta$ are hardcoded in many cases.**
- folder ./rcpp/ contains `rcpp` source files for posteriors (up to a constant) and gradients
- folder ./sims/ contains R code to reproduce all simulations and plots from the paper and the supplement
- folder ./framingham/ contains R code to reproduce all models and plots for the Framingham heart study illustration
- folder ./leukemia/ contains R code to reproduce all models and plots for the Leukemia lung cancer survival illustration
- folder ./veteran/ contains R code to reproduce all models and plots for the Veteran illustration

Used libraries are captured in the renv.lock file (after installing the `renv`package use `renv::init()`).

### Usage
All code in sims/framingham/leukemia/veteran is self-contained and can be used by simply running the scripts.


## Processed Files
Folder ./processed_data/ contains all results used in the paper except for files that are larger than GitHub's limits.

**Importantly, the authors should provide an overview of how to carry
out the analyses presented in their manuscript in the `README.md` of their
repository, replacing the content in this file.** This overview would
generally refer to scripts/code files that execute the analyses and are
placed either in the main directory or the `/code` subdirectory. The
*Workflow* section of the ACC form should refer to this README.md as
containing the instructions for how to reproduce the analyses.

### Step 3

Author(s) use `git commit` to track changes over time and use `git push`
to push changes to a repository on the author(s) personal GitHub
account.

### Step 4

Author(s) submit a link to their GitHub repository as part of the [JASA
Reproducibility review process](https://jasa-acs.github.io/repro-guide/),
required upon submission of an invited revision.

### Step 5

JASA Associate Editors for Reproducibility will review the materials in
the GitHub repository of the authors and submit a
reproducibility review as part of the standard JASA review process.
Authors have the opportunity to respond to the review by making changes
and pushing their changes to their personal GitHub repository.

### Step 6

Once the manuscript is accepted, the materials in the author(s) personal
GitHub repository will be copied to the [JASA repository](https://github.com/jasa-acs).

## Reproducibility materials file structure

This template provides a suggested file structure for a JASA submission, but authors are free
to modify this structure.

The suggested components are as follows. Directories in the submission may have subdirectories to
further organize the materials.

1.  A `README.md` file - This file gives a short description of the
    paper and an overview of how to carry out the analyses presented in their manuscript.
2.  A `manuscript` directory - This directory will generally hold the source files
    (often LaTeX or Rmd) for the manuscript and any files directly related to the
    generation of the manuscript, including figure files.
3.  A `data` directory - This directory will generally hold the real data files 
    (or facsimile versions of them in place of confidential data) and simulated data files.
    See `data/README.md` for more details. 
4.  A `code` directory - This directory will generally hold 
    source code files that contain the core code to implement the method and various utility/auxiliary functions.
5.  An `output` directory - This directory will generally hold objects derived
    from computations, including results of simulations or real data analyses. See `output/README.md` for more details.

## Guidance on the use of reproducible environments

Submissions may include the use of reproducible environments capturing
state of a machine generating manuscript artifacts and even the
manuscript itself. Here we discuss two types of reproducible
environments and their use. Both virtual and package environments may be
put in the `code` directory.

### Package environments

Package environments capture the set of packages used by a programming
language needed to generate output. The R programming language has
`renv`, `switchr` and others to accomplish this, Python has `venv`,
`conda` and others, and Julia has native support (through the `Pkg`
package). When submitting these types of environments, the following are
suggested.

1.  Clearly indicate (in the overall `README.md`) the language(s) used (including version) 
    and the package environment tool used (e.g., `renv`, `conda`).
2.  Use a single package environment for all reproducible content.
3.  Prefer packages from package archives (CRAN, Bioconductor,
    RForge.net for example).
4.  If you use packages from a code repository (GitHub, GitLab, etc.)
    then use a release version if possible, or indicate the commit used. You could also consider
    forking the repository and providing a release.

### Virtual environments

Virtual environments such as Docker and Singlarity capture
the entire computing environment in which computations were performed.
In general, they are a more robust solution, capable of taking a
“snapshot” of a machine, including any system-level utilities and
external libraries needed to perform your computations. They have the
advantage that reproducing materials means running the virtual
environment, rather than recreating the programming language environment.
If using a virtual environment, we ask that 
you provide a definition file (e.g., a Dockerfile) or (perhaps better)
a link to an image in a standard online registry, such as DockerHub.

## References

Gentleman, Robert, and Duncan Temple Lang. “[Statistical Analyses and
Reproducible
Research](http://biostats.bepress.com/cgi/viewcontent.cgi?article=1001&context=bioconductor).”
(2004).

Gentleman, Robert. “[Reproducible research: a bioinformatics case
study](https://www.degruyter.com/document/doi/10.2202/1544-6115.1034/html).”
Statistical applications in genetics and molecular biology 4.1 (2005).

Marwick, Ben, and Bryan, Jennifer, and Attali, Dean, and Hollister,
Jeffrey W. [rrrpkg Github Page](https://github.com/ropensci/rrrpkg).
