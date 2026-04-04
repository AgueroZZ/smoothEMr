# Collapsed-ML csmoothEM with weighted E-step and weighted hyper-step

Collapsed-ML csmoothEM with weighted E-step and weighted hyper-step

## Usage

``` r
do_csmoothEM_ml_collapsed_weighted(
  object,
  data,
  feature_weights,
  iter = 1,
  record = TRUE,
  lambda_min = NULL,
  lambda_max = NULL,
  sigma_update = c("ml", "mstep"),
  verbose = FALSE
)
```
