# Evaluate soft partition accuracy (handles label-switching)

Evaluate soft partition accuracy (handles label-switching)

## Usage

``` r
evaluate_partition(result, true_assign)
```

## Arguments

- result:

  Output of soft_two_trajectory_EM().

- true_assign:

  Ground-truth vector ("A","B","noise").

## Value

List: accuracy, best_alignment, confusion_table.
