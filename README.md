
# IVPP

<!-- badges: start -->
[![R-CMD-check](https://github.com/xinkaidupsy/IVPP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/xinkaidupsy/IVPP/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

The goal of IVPP is to compare network models for intensive time-series and panel data. 

## Installation

You can install the development version of IVPP from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xinkaidupsy/IVPP")
```

## Example

An example that uses IVPP to compare panelGVAR models:

``` r
library(IVPP)
# Generate the network
net_ls <- gen_panelGVAR(n_node = 6,
                        p_rewire = 0.5,
                        n_group = 3)

# Generate the data
data <- sim_panelGVAR(temp_base_ls = net_ls$temporal,
                      cont_base_ls = net_ls$omega_zeta_within,
                      n_person = 500,
                      n_time = 4,
                      n_group = 3,
                      n_node = 6)

# IVPP test on the temporal network
ivpp <- IVPP_panelgvar(data,
                       vars = paste0("V",1:6),
                       idvar = "subject",
                       beepvar = "time",
                       groups = "group",
                       test = "temporal",
                       net_type = "saturated",
                       prune_net = "temporal",
                       partial_prune = TRUE,
                       estimator = "FIML",
                       standardize = "z")

```

An example that uses IVPP to compare N = 1 GVAR models

``` r
library(IVPP)

# Generate the network
net_ls <- gen_tsGVAR(n_node = 6,
                     p_rewire = 0.5,
                     n_persons = 3)

# Generate the data
data <- sim_tsGVAR(beta_base_ls = net_ls$beta,
                   kappa_base_ls = net_ls$kappa,
                   # n_person = 3,
                   n_time = 50)

# IVPP test on
ivpp_ts <- IVPP_tsgvar(data = data,
                       vars = paste0("V",1:6),
                       idvar = "id",
                       test = "temporal",
                       net_type = "saturated",
                       prune_net = "temporal",
                       partial_prune = TRUE,
                       estimator = "FIML",
                       standardize = "z")
```
