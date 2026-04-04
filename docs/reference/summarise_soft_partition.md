# Summarise a soft partition result

Summarise a soft partition result

## Usage

``` r
summarise_soft_partition(result, feature_names = NULL)
```

## Arguments

- result:

  Output of soft_two_trajectory_EM().

- feature_names:

  Optional character vector length d.

## Value

data.frame with feature, w_A, w_B, assign, entropy.
