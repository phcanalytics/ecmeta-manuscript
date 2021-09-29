# Benchmark hazard ratios ------------------------------------------------------
#' Benchmark hazard ratios
#' 
#' Obtain estimated hazard ratios in a tidy data frame from a `surv_ate` object and
#' compare them to benchmark hazard ratios (e.g., from target randomized clinical
#' trials). Quickly plot comparisons of the estimated hazard ratios to the benchmarks
#' using `autoplot()`. Generate a pretty HTML table with 95 percent confidence 
#' intervals using `html_table()`.
#' 
#' @param object A [`surv_ate`] object.
#' @param benchmarks A data frame containing benchmark hazard ratios, potentially
#' by group. It must 
#' @param becnhmark_estimate The column identifier for the benchmark hazard ratio
#' estimate.
#' @param benchmark_lower The column identifier for the lower confidence limit
#' of the benchmark hazard ratio.
#' @param benchmark_upper The column identifier for the upper confidence limit
#' of the benchmark hazard ratio.
#' 
#' @return A `hazard_ratio` object that inherits from [`dplyr::tibble`].
#' @export
benchmark_hazard_ratios <- function(object, benchmarks = NULL, benchmark_estimate = "estimate",
                                    benchmark_lower = "lower", benchmark_upper = "upper") {
  
  # Tidy hazard ratios
  hr <- hazard_ratio(object)
  
  # Add benchmark estimates
  b <- copy(benchmarks)
  data.table::setnames(b, benchmark_estimate, "benchmark_estimate")
  data.table::setnames(b, benchmark_lower, "benchmark_lower")
  data.table::setnames(b, benchmark_upper, "benchmark_upper")
  if (inherits(object, "grouped_df")) {
    hr <-  dplyr::left_join(hr, b, by = dplyr::group_vars(object))
  } else {
    hr <- cbind(hr, b)
  }
  
  
  # Return
  class(hr) <- c("benchmarked_hazard_ratios", class(hr))
  attr(hr, "treat") <- attr(object, "treat")
  hr
}

#' @rdname benchmark_hazard_ratios
#' @export
autoplot.benchmarked_hazard_ratios <- function(object,
                                               xlab = "Benchmark hazard ratio",
                                               ylab = "Estimated hazard ratio",
                                               scale_breaks = seq(.4, 2, .2),
                                               ...) {
  
  object$method_label <- label_ps_method(object$method)
  object <- object[object$term == attr(object, "treat"), ]
  
  # Axis limits
  est_min <- min(object[, c("estimate", "benchmark_estimate")])
  est_max <- max(object[, c("estimate", "benchmark_estimate")])
  lim_lower <-  floor(est_min * 10) /10
  lim_upper <- ceiling(est_max * 10)/10
  lims <- c(lim_lower, lim_upper)
  if (is.null(scale_breaks)) scale_breaks <- seq(lim_lower, lim_upper, 
                                                 length.out = 20)
  
  # Plot
  grp_name <- dplyr::group_vars(object)
  ggplot(object, aes(x = benchmark_estimate, y = estimate)) +
    geom_point() +
    facet_wrap(~method_label) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    xlab(xlab) +
    ylab(ylab) +
    scale_x_continuous(breaks = scale_breaks, limits = lims) +
    scale_y_continuous(breaks = scale_breaks, limits = lims) +
    geom_text(aes(label = .data[[grp_name]]), hjust = 0, nudge_x = 0.015,
              size = 3)
}

#' @rdname benchmark_hazard_ratios
#' @export
html_table.benchmarked_hazard_ratios <- function(x, ...) {
  treat <- attr(x, "treat")
  grp_name <- dplyr::group_vars(x)
  x2 <- as.data.table(x)
  x2 <- x2[term == treat]
  
  # Format confidence intervals
  add_ci <- function(x, l, u, digits = 2) {
    paste0(
      formatC(x, format = "f", digits = digits),
      " (",
      formatC(l, format = "f", digits = digits),
      "-",
      formatC(u, format = "f", digits = digits),
      ")"
    )
  }
  
  x2[, festimate := add_ci(estimate, lower, upper)]
  x2[, fbenchmark := add_ci(benchmark_estimate, benchmark_lower, benchmark_upper)]
  
  # Convert to wide format
  f_dcast <- as.formula(paste0(grp_name, " + term + fbenchmark ~", 
                               " method"))
  
  x2 <- dcast(as.data.table(x2), formula = f_dcast, value.var = "festimate")
  data.table::setnames(x2, "fbenchmark", "benchmark")
  
  # Make HTML table
  knitr::kable(x2) %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "responsive"))
}