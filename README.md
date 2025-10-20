
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dist.pois.nbinom.hamim

<!-- badges: start -->

<!-- badges: end -->

dist.pois.nbinom.hamim provides tools to preview and test chi-square
goodness-of-fit (GOF) for Poisson and Negative Binomial on discrete
count data, with two-sided automatic pooling (so each final class has Ex
≥ threshold) and an optional single tail bin (“t+”). Multiple estimation
methods for NegBin are supported (MLE, direct-k, manual, log-iterative).

## Installation

# From a local checkout (recommended during development)

# Run in the package root (where DESCRIPTION lives)

devtools::install(build_vignettes = TRUE)

# Or, once you publish to GitHub, replace USER/REPO accordingly:

# remotes::install_github(“USER/dist.pois.nbinom.hamim”, build_vignettes = TRUE)

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

## Example

library(dist.pois.nbinom.hamim)

# Example data

txt \<- ” X Fx 0 114 1 25 2 15 3 10 4 6 5 5 6 2 7 1 8 1 9 0 10 1 ” dat
\<- read.table(text = txt, header = TRUE)

# 1) Poisson preview (no pooling), include tail “t+”

prev \<- distribusi.hamim( data = dat, distrib = “poisson”, pool =
“none”, tail = TRUE, plot = FALSE ) prev$tabel
sum(prev$tabel\$Px) \# ~1 because tail bin is included

# 2) Poisson GOF with automatic pooling (Ex \>= 3)

gof_pois \<- distribusi.hamim( data = dat, distrib = “poisson”, ex_min =
3, pool = “auto”, tail = TRUE, plot = FALSE )
gof_pois$p_val; gof_pois$df head(gof_pois\$tabel)

# 3) Negative Binomial (log.iterasi), Ex \>= 5

gof_nb \<- distribusi.hamim( data = dat, distrib = “negbinom”, method =
“log.iterasi”, ex_min = 5, pool = “auto”, tail = TRUE, plot = FALSE )
gof_nb$k; gof_nb$mu; gof_nb\$p_val

``` r
library(dist.pois.nbinom.hamim)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.

info.distribusi.hamim() \# short reminders (ID, GI, d-stat, GOF) \#
browseVignettes(“dist.pois.nbinom.hamim”) \# open the vignette
