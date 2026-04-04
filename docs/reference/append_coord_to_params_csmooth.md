# Append a coordinate to csmooth-style parameters (internal)

Appends a single coordinate (feature) to an existing csmooth-style
parameter list. This is used in greedy feature assignment to grow a
partition without reinitializing from scratch.

## Usage

``` r
append_coord_to_params_csmooth(params, one)
```

## Arguments

- params:

  List with fields `pi`, `mu`, `sigma2` for the current partition. `mu`
  must be a list of length K, each element a numeric vector of length
  d_sub. `sigma2` must be a numeric vector of length d_sub
  (homoskedastic case).

- one:

  List describing the new coordinate. Must contain:

  - `mu_vec`: numeric vector of length K (component means for the new
    coordinate)

  - `sigma2`: numeric scalar (variance for the new coordinate)

## Value

Updated `params` with the new coordinate appended.
