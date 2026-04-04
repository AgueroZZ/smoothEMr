# Drop one coordinate from csmooth-style parameters (internal)

Removes a single coordinate (feature) from a csmooth-style parameter
list. Used for warm-start backward greedy moves (remove from one
partition).

## Usage

``` r
drop_coord_from_params_csmooth(params, pos)
```

## Arguments

- params:

  List with fields `pi`, `mu`, `sigma2`. `mu` must be a list of length
  K, each element length d_sub. `sigma2` must be a numeric vector length
  d_sub (homoskedastic).

- pos:

  Integer in `1:d_sub`, position of the coordinate to remove (in the
  CURRENT partition's feature space).

## Value

Updated `params`.
