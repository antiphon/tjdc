# `tjdc`: R package for Trends and Jumps in (flat) DataCubes


## Purpose

This R-package accompanies the research paper

> T Rajala, P Packalen, M Myllymäki, A Kangas (2023): Improving detection of changepoints in short and noisy time-series with local correlations: Connecting the events in pixel neighbourhoods, *Journal of Agricultural, Biological and Environmental Statistics*, https://doi.org/10.1007/s13253-023-00546-1
 
and implements the (proof-of-concept) Gibbs sampler for the described model. 

The described model was developed because no ready made solution was found for the required analysis; see below.


## Installation

In addition to CRAN dependencies, install the `looptimer` package.

```r
# 
devtools::install_github("antiphon/looptimer")

# install.packages("doParallel", "foreach", "sf", "stars", "dplyr", "tidyr"))
devtools::install_github("antiphon/tjdc")
```

This package has a vignette, so consider using `build_vignettes = TRUE`. The vignette is also readable at the package site, 


## Input data

Consider a timeseries data, stored as a vector $y_i=(y_{i1},...,y_{iT})$ of values 
$y_{it}\in\mathbf{R}$ recorded at times $t_1 < ... < t_T$. 

Now consider that each timeseries $y_i$ is recorded on the nodes of (spatial) grid. Such as grid-time data array can be called a *datacube*: Dimension x-y of the cuboid describe location, and z describes time. 

See the example datasets

```r
data(package = "tjdc")

data("case_medium_zoom1")
data("case_medium_zoom1_details")
data("easystripes")
```

The datasets analysed in the original paper are available at Zenodo, https://zenodo.org/record/8009800 .

## The model

The package, or more specifically the model it implements, is designed to analyse datacubes assuming

1. We want to know if, and when, there is a radical change in the mean of each $y_i$
2. The series are short, $T\ll 100$, like $T=11$ in the paper's main example
3. The unstructed variability (noise) is relatively high, e.g. $CV \gg 5\%$
4. The changepoint events are spatially correlated, so taking place in "patches" of locations

## Functionality

The package provides wrappers for 

1. Running the Gibbs sampler
1. Parallelising: Splitting the datacube spatially with minor overlap, parallel running the model, and then gathering the posteriors

The main fitting function is `tjdc::tj_fit_m0.6`, proceed there.

For large problems it is recommended to use the file-storage options (see `tjdc::tj_fit_m0.6`).



