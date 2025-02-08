
# IVPP

<p align="center">
    <a href="https://www.repostatus.org/#active"><img src="https://www.repostatus.org/badges/latest/active.svg" alt="Repository status"/></a>
    <a href="https://github.com/mihaiconstantin/powerly/releases"><img src="https://img.shields.io/github/v/release/mihaiconstantin/powerly?display_name=tag&sort=semver"/></a>
    <a href="https://www.r-pkg.org/pkg/powerly"><img src="https://www.r-pkg.org/badges/version/powerly" alt="CRAN version"/></a>
    <a href="https://www.r-pkg.org/pkg/powerly"><img src="https://cranlogs.r-pkg.org/badges/grand-total/powerly" alt="CRAN RStudio mirror downloads"/></a>
    <a href="https://app.codecov.io/gh/mihaiconstantin/powerly"><img src="https://codecov.io/gh/mihaiconstantin/powerly/branch/main/graph/badge.svg?token=YUCO8ULBCM" alt="Code coverage"/></a>
    <a href="https://github.com/mihaiconstantin/powerly/actions"><img src="https://github.com/mihaiconstantin/powerly/workflows/R-CMD-check/badge.svg" alt="R-CMD-check" /></a>
    <a href="https://cran.rstudio.com/web/checks/check_results_powerly.html"><img src="https://badges.cranchecks.info/worst/powerly.svg" alt="CRAN checks"/></a>
</p>


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
