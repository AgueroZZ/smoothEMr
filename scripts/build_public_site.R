#!/usr/bin/env Rscript

override <- list(home = list(sidebar = FALSE))

public_reference_topics <- c(
  "MPCurver",
  "mpcurve",
  "fit_mpcurve",
  "do_mpcurve",
  "fitted_prior",
  "print.mpcurve",
  "summary.mpcurve",
  "plot.mpcurve",
  "PCA_ordering",
  "fiedler_ordering",
  "isomap_ordering",
  "pcurve_ordering",
  "tSNE_ordering",
  "simulate_cavi_toy",
  "simulate_dual_trajectory",
  "simulate_intrinsic_trajectories",
  "simulate_spiral2d",
  "simulate_swiss_roll_1d_2d",
  "simulate_two_order_gp_dataset"
)

public_articles <- c(
  "mpcurve_intro",
  "partition",
  "fitness"
)

rewrite_pkgdown_metadata <- function(path, public_articles) {
  if (!file.exists(path)) {
    return(invisible(NULL))
  }

  lines <- readLines(path, warn = FALSE)
  output <- character()
  in_articles <- FALSE

  for (line in lines) {
    if (identical(line, "articles:")) {
      in_articles <- TRUE
      output <- c(output, line)
      next
    }

    if (in_articles && grepl("^[^[:space:]]", line)) {
      in_articles <- FALSE
    }

    if (in_articles) {
      match <- regmatches(
        line,
        regexec("^  ([^:]+): [^[:space:]]+\\.html$", line)
      )[[1]]

      if (length(match) > 1L) {
        if (match[2] %in% public_articles) {
          output <- c(output, line)
        }
        next
      }
    }

    output <- c(output, line)
  }

  writeLines(output, path)
}

if (dir.exists("docs")) {
  unlink("docs", recursive = TRUE, force = TRUE)
}

dir.create("docs", recursive = TRUE, showWarnings = FALSE)

pkg_home <- pkgdown:::section_init(".", override = override)
pkgdown:::build_home_index(pkg_home, quiet = TRUE)
if (!pkg_home$development$in_dev) {
  pkgdown:::build_404(pkg_home)
}

pkg_ref <- pkgdown:::section_init(".", "reference", override = override)
pkg_ref$topics <- pkg_ref$topics[
  pkg_ref$topics$name %in% public_reference_topics,
  ,
  drop = FALSE
]
pkgdown:::build_reference_index(pkg_ref)
pkgdown::build_reference(
  pkg = pkg_ref,
  lazy = FALSE,
  examples = FALSE,
  preview = FALSE,
  devel = FALSE
)

pkg_art <- pkgdown:::section_init(".", "articles", override = override)
pkg_art$vignettes <- pkg_art$vignettes[
  pkg_art$vignettes$name %in% public_articles,
  ,
  drop = FALSE
]
pkgdown:::build_articles_index(pkg_art)
for (name in public_articles) {
  pkgdown:::build_article(
    name,
    pkg = pkg_art,
    lazy = FALSE,
    seed = 1014L,
    new_process = FALSE,
    quiet = TRUE
  )
}

pkgdown:::build_search(pkg_home)
pkgdown:::build_sitemap(pkg_home)
rewrite_pkgdown_metadata("docs/pkgdown.yml", public_articles)
