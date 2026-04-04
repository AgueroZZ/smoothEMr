# Package index

## All functions

- [`EM_algorithm()`](https://aguerozz.github.io/MPCurver/reference/EM_algorithm.md)
  : Penalized Expectation-Maximization (EM) algorithm for SmoothEM

- [`PCA_ordering()`](https://aguerozz.github.io/MPCurver/reference/PCA_ordering.md)
  : PCA ordering (PC1)

- [`as_csmooth_em()`](https://aguerozz.github.io/MPCurver/reference/as_csmooth_em.md)
  : Construct a csmooth_em object

- [`as_mpcurve()`](https://aguerozz.github.io/MPCurver/reference/as_mpcurve.md)
  :

  Convert an algorithm fit object to an `mpcurve` model object

- [`as_smooth_em()`](https://aguerozz.github.io/MPCurver/reference/as_smooth_em.md)
  : Construct a smooth_em object from an EM_algorithm fit

- [`backward_two_ordering_partition_csmooth()`](https://aguerozz.github.io/MPCurver/reference/backward_two_ordering_partition_csmooth.md)
  : Backward greedy feature partition into two orderings (csmoothEM,
  warm-start only)

- [`cache_csmooth_params()`](https://aguerozz.github.io/MPCurver/reference/cache_csmooth_params.md)
  : Cache inverse-sigma2 and log-determinant for csmooth_em (diagonal
  covariance)

- [`cavi()`](https://aguerozz.github.io/MPCurver/reference/cavi.md) :
  Fit the recommended CAVI model for a single MPCurver ordering

- [`compare_positions_by_history()`](https://aguerozz.github.io/MPCurver/reference/compare_positions_by_history.md)
  : Compare posterior position summaries across two histories

- [`do_cavi()`](https://aguerozz.github.io/MPCurver/reference/do_cavi.md)
  :

  Continue CAVI sweeps on an existing `cavi` fit

- [`do_csmoothEM()`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM.md)
  : Run csmoothEM iterations

- [`do_csmoothEM_ml_collapsed()`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM_ml_collapsed.md)
  : Run collapsed-ML csmoothEM iterations (homoskedastic only)

- [`do_csmoothEM_weighted()`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM_weighted.md)
  : Run csmoothEM with weighted inner updates

- [`do_mpcurve()`](https://aguerozz.github.io/MPCurver/reference/do_mpcurve.md)
  :

  Continue CAVI iterations on an existing `mpcurve` fit

- [`do_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/do_smoothEM.md)
  : Run SmoothEM for a given number of iterations on a smooth_em object

- [`evaluate_partition()`](https://aguerozz.github.io/MPCurver/reference/evaluate_partition.md)
  : Evaluate soft partition accuracy (handles label-switching)

- [`fiedler_ordering()`](https://aguerozz.github.io/MPCurver/reference/fiedler_ordering.md)
  : Fiedler ordering from a kNN graph

- [`fit_mpcurve()`](https://aguerozz.github.io/MPCurver/reference/fit_mpcurve.md)
  : Fit an MPCurve model with the public CAVI interface

- [`forward_two_ordering_partition_csmooth()`](https://aguerozz.github.io/MPCurver/reference/forward_two_ordering_partition_csmooth.md)
  : Forward greedy feature partition into two orderings (csmoothEM
  version)

- [`generalized_logdet()`](https://aguerozz.github.io/MPCurver/reference/generalized_logdet.md)
  : Generalized log-determinant (sum of log positive eigenvalues)

- [`get_history_block()`](https://aguerozz.github.io/MPCurver/reference/get_history_block.md)
  : Extract a (u, gamma) block for a given progressive history index

- [`greedy_backward_filter_csmooth()`](https://aguerozz.github.io/MPCurver/reference/greedy_backward_filter_csmooth.md)
  : Greedy backward filtering of features for csmoothEM ordering
  inference

- [`init_cov_cache_fast()`](https://aguerozz.github.io/MPCurver/reference/init_cov_cache_fast.md)
  : Cache covariance inverses and log-determinants

- [`init_m_trajectories_cavi()`](https://aguerozz.github.io/MPCurver/reference/init_m_trajectories_cavi.md)
  : Initialize M CAVI trajectories for multi-ordering partitioning

- [`init_two_trajectories()`](https://aguerozz.github.io/MPCurver/reference/init_two_trajectories.md)
  : Initialise two trajectory fits guaranteed to face different
  directions

- [`init_two_trajectories_cavi()`](https://aguerozz.github.io/MPCurver/reference/init_two_trajectories_cavi.md)
  : Initialize two CAVI trajectories for dual-ordering partitioning

- [`initialize_csmoothEM()`](https://aguerozz.github.io/MPCurver/reference/initialize_csmoothEM.md)
  : Initialize csmoothEM (coordinate-specific SmoothEM)

- [`initialize_ordering()`](https://aguerozz.github.io/MPCurver/reference/initialize_ordering.md)
  : Initialize GMM parameters using an ordering method

- [`initialize_ordering_csmooth()`](https://aguerozz.github.io/MPCurver/reference/initialize_ordering_csmooth.md)
  : Initialize ordering for csmoothEM (diagonal-variance version)

- [`initialize_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/initialize_smoothEM.md)
  : Initialize SmoothEM (single- or multi-scale)

- [`isomap_ordering()`](https://aguerozz.github.io/MPCurver/reference/isomap_ordering.md)
  : Landmark Isomap ordering (no vegan; avoids O(n^2) dist matrix)

- [`krige_mu_list_to_full_grid()`](https://aguerozz.github.io/MPCurver/reference/krige_mu_list_to_full_grid.md)
  : Krige SmoothEM mean functions from an observed grid to the final
  grid

- [`kriging_from_precision()`](https://aguerozz.github.io/MPCurver/reference/kriging_from_precision.md)
  : Krige/interpolate from observed nodes using a Gaussian precision
  matrix

- [`lambda_scale_for_spacing()`](https://aguerozz.github.io/MPCurver/reference/lambda_scale_for_spacing.md)
  : Scale random-walk penalty strength to account for grid spacing

- [`logsumexp()`](https://aguerozz.github.io/MPCurver/reference/logsumexp.md)
  : Numerically stable log-sum-exp

- [`make_VAR1_precision()`](https://aguerozz.github.io/MPCurver/reference/make_VAR1_precision.md)
  : Prior precision for VAR(1) on K time points (d-dimensional)

- [`make_default_init()`](https://aguerozz.github.io/MPCurver/reference/make_default_init.md)
  : Default random initialization for a Gaussian mixture model

- [`make_hierarchical_levels()`](https://aguerozz.github.io/MPCurver/reference/make_hierarchical_levels.md)
  : Construct a hierarchy of nested grids inside a final grid

- [`make_init()`](https://aguerozz.github.io/MPCurver/reference/make_init.md)
  : Initialization from an ordering vector

- [`make_init_csmooth()`](https://aguerozz.github.io/MPCurver/reference/make_init_csmooth.md)
  : Initialization from an ordering vector for csmoothEM

- [`make_lattice_rwq_precision()`](https://aguerozz.github.io/MPCurver/reference/make_lattice_rwq_precision.md)
  : Prior precision for q-th order random walk on a K x K lattice
  (d-dimensional)

- [`make_lattice_rwq_precision_sparse()`](https://aguerozz.github.io/MPCurver/reference/make_lattice_rwq_precision_sparse.md)
  : Sparse prior precision for q-th order random walk on a K x K lattice
  (d-dimensional)

- [`make_random_walk_precision()`](https://aguerozz.github.io/MPCurver/reference/make_random_walk_precision.md)
  : Prior precision for q-th order random walk on K time points
  (d-dimensional)

- [`make_random_walk_precision_sparse()`](https://aguerozz.github.io/MPCurver/reference/make_random_walk_precision_sparse.md)
  : Sparse prior precision for q-th order random walk on K time points
  (d-dimensional)

- [`match_locations_to_grid()`](https://aguerozz.github.io/MPCurver/reference/match_locations_to_grid.md)
  : Match locations on an old grid to indices on a new grid

- [`optimize_initial()`](https://aguerozz.github.io/MPCurver/reference/optimize_initial.md)
  : Choose the best SmoothEM initialization by ELBO

- [`optimize_initial_csmoothEM()`](https://aguerozz.github.io/MPCurver/reference/optimize_initial_csmoothEM.md)
  : Optimize csmoothEM initialization by comparing multiple methods
  after a warm start

- [`parallel_initial()`](https://aguerozz.github.io/MPCurver/reference/parallel_initial.md)
  : Run multiple SmoothEM initializations in parallel

- [`parallel_initial_csmoothEM()`](https://aguerozz.github.io/MPCurver/reference/parallel_initial_csmoothEM.md)
  : Run multiple csmoothEM initializations in parallel and summarize
  results

- [`partition_features_twofits()`](https://aguerozz.github.io/MPCurver/reference/partition_features_twofits.md)
  : Partition features by comparing per-coordinate collapsed scores from
  two fits

- [`partition_features_twofits_cavi()`](https://aguerozz.github.io/MPCurver/reference/partition_features_twofits_cavi.md)
  : Partition features by comparing CAVI feature scores from two fits

- [`pcurve_ordering()`](https://aguerozz.github.io/MPCurver/reference/pcurve_ordering.md)
  : Principal curve ordering

- [`plot(`*`<cavi>`*`)`](https://aguerozz.github.io/MPCurver/reference/plot.cavi.md)
  :

  Plot a `cavi` fit

- [`plot(`*`<csmooth_em>`*`)`](https://aguerozz.github.io/MPCurver/reference/plot.csmooth_em.md)
  : Plot a csmooth_em object

- [`plot(`*`<mpcurve>`*`)`](https://aguerozz.github.io/MPCurver/reference/plot.mpcurve.md)
  :

  Plot an `mpcurve` model fit

- [`plot(`*`<smooth_em>`*`)`](https://aguerozz.github.io/MPCurver/reference/plot.smooth_em.md)
  : Plot a smooth_em object

- [`plot_EM_embedding()`](https://aguerozz.github.io/MPCurver/reference/plot_EM_embedding.md)
  : Plot data embedding with SmoothEM component means

- [`plot_EM_embedding2D()`](https://aguerozz.github.io/MPCurver/reference/plot_EM_embedding2D.md)
  : Plot a 2D embedding with Smooth-EM component means overlaid

- [`plot_coordinate_change()`](https://aguerozz.github.io/MPCurver/reference/plot_coordinate_change.md)
  : Plot change in a single coordinate between two histories

- [`plot_mu_history()`](https://aguerozz.github.io/MPCurver/reference/plot_mu_history.md)
  : Plot kriged mean curves from mu_full_history

- [`plot_order_EM_overlay2D()`](https://aguerozz.github.io/MPCurver/reference/plot_order_EM_overlay2D.md)
  : Plot a 2D embedding colored by an ordering vector with EM means
  overlaid

- [`plot_soft_weights()`](https://aguerozz.github.io/MPCurver/reference/plot_soft_weights.md)
  : Plot soft feature weights (4-panel diagnostic)

- [`print(`*`<csmooth_em>`*`)`](https://aguerozz.github.io/MPCurver/reference/print.csmooth_em.md)
  : Print csmooth_em object

- [`print(`*`<smooth_em>`*`)`](https://aguerozz.github.io/MPCurver/reference/print.smooth_em.md)
  : Print method for smooth_em objects

- [`print(`*`<summary.csmooth_em>`*`)`](https://aguerozz.github.io/MPCurver/reference/print.summary.csmooth_em.md)
  : Print summary of csmooth_em

- [`print(`*`<summary.smooth_em>`*`)`](https://aguerozz.github.io/MPCurver/reference/print.summary.smooth_em.md)
  : Print method for summary.smooth_em objects

- [`progressive_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/progressive_smoothEM.md)
  : Progressive-resolution initialization for SmoothEM via kriging on a
  final grid

- [`score_feature_given_Gamma()`](https://aguerozz.github.io/MPCurver/reference/score_feature_given_Gamma.md)
  : Score a single feature given responsibilities (Gamma)

- [`score_features_onefit()`](https://aguerozz.github.io/MPCurver/reference/score_features_onefit.md)
  : Score features under a fitted csmooth_em model via per-coordinate
  collapsed contributions

- [`score_features_onefit_cavi()`](https://aguerozz.github.io/MPCurver/reference/score_features_onefit_cavi.md)
  :

  Score features under a fitted `cavi` model

- [`simulate_cavi_toy()`](https://aguerozz.github.io/MPCurver/reference/simulate_cavi_toy.md)
  : Simulate a toy dataset from the current single-ordering CAVI model

- [`simulate_dual_trajectory()`](https://aguerozz.github.io/MPCurver/reference/simulate_dual_trajectory.md)
  : Simulate a dual-trajectory dataset with configurable trajectory
  families

- [`simulate_intrinsic_trajectories()`](https://aguerozz.github.io/MPCurver/reference/simulate_intrinsic_trajectories.md)
  : Simulate a multi-ordering intrinsic-trajectory dataset

- [`simulate_spiral2d()`](https://aguerozz.github.io/MPCurver/reference/simulate_spiral2d.md)
  : Simulate a 2D Archimedean spiral (helix-like) with noise

- [`simulate_swiss_roll_1d_2d()`](https://aguerozz.github.io/MPCurver/reference/simulate_swiss_roll_1d_2d.md)
  : Simulate a 2D "swiss roll" spiral with a 1D latent parameter

- [`simulate_two_order_gp_dataset()`](https://aguerozz.github.io/MPCurver/reference/simulate_two_order_gp_dataset.md)
  : Simulate a two-ordering GP dataset (Matern) for feature partitioning

- [`soft_partition_cavi()`](https://aguerozz.github.io/MPCurver/reference/soft_partition_cavi.md)
  : Soft M-ordering partition using mean-field CAVI

- [`soft_two_trajectory_EM()`](https://aguerozz.github.io/MPCurver/reference/soft_two_trajectory_EM.md)
  : Soft dual-trajectory partition via annealed CAVI

- [`soft_two_trajectory_cavi()`](https://aguerozz.github.io/MPCurver/reference/soft_two_trajectory_cavi.md)
  : Soft dual-trajectory partition using explicit mean-field CAVI

- [`subset_csmooth_em_fit()`](https://aguerozz.github.io/MPCurver/reference/subset_csmooth_em_fit.md)
  : Subset a csmooth_em object by features (columns)

- [`summarise_soft_partition()`](https://aguerozz.github.io/MPCurver/reference/summarise_soft_partition.md)
  : Summarise a soft partition result

- [`summary(`*`<cavi>`*`)`](https://aguerozz.github.io/MPCurver/reference/summary.cavi.md)
  :

  Summary method for `cavi` objects

- [`summary(`*`<csmooth_em>`*`)`](https://aguerozz.github.io/MPCurver/reference/summary.csmooth_em.md)
  : Summary for csmooth_em object

- [`summary(`*`<mpcurve>`*`)`](https://aguerozz.github.io/MPCurver/reference/summary.mpcurve.md)
  :

  Summarise an `mpcurve` model fit

- [`summary(`*`<smooth_em>`*`)`](https://aguerozz.github.io/MPCurver/reference/summary.smooth_em.md)
  : Summary method for smooth_em objects

- [`tSNE_ordering()`](https://aguerozz.github.io/MPCurver/reference/tSNE_ordering.md)
  : t-SNE ordering
