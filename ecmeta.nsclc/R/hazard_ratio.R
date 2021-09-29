#' Hazard ratios
#' 
#' Obtain estimated hazard ratios in a tidy data frame from a `surv_ate` object and
#' quickly plot them using `autoplot()`.
#' 
#' @param object A [`surv_ate`] object.
#' @param exponentiate Logical indicating whether or not to exponentiate
#' the coefficient estimates. Default is `TRUE` meaning that hazard ratios
#' are returned; if `FALSE`, then log hazard ratios are returned.
#' 
#' @return A `hazard_ratio` object that inherits from [`dplyr::tibble`].
#' @export
hazard_ratio <- function(object, exponentiate = TRUE) {
  if (!inherits(object, c("surv_ate", "grouped_surv_ate"))) {
    stop("'object' must be of class 'surv_ate' or 'grouped_surv_ate'.",
         call. = FALSE)
  }
  
  loghr_list <- object$loghr
  names(loghr_list) <- object$method
  hr <- dplyr::bind_rows(loghr_list, .id = "method")
  if (exponentiate) {
    for (v in c("estimate", "lower", "upper")) {
      hr[, v] <- exp(hr[, v])
    }
  }

  # Add group identifier if required
  if (inherits(object, "grouped_surv_ate")) {
    grp_name <- dplyr::group_vars(object)
    n_rep <- sapply(object$loghr, nrow)
    grp_df <- object[, grp_name]
    hr <- cbind(
      grp_df[rep(seq_len(nrow(grp_df)), times = n_rep), ],
      hr
    )
  }
  
  # Limit to treatment
  hr <- hr[hr$term == attr(object, "treat"), ]
  
  # Return
  class(hr) <- c("hazard_ratio", class(hr))
  attr(hr, "treat") <- attr(object, "treat")
  attr(hr, "exponentiate") <- exponentiate
  hr
}

#' @rdname hazard_ratio
#' @export
autoplot.hazard_ratio <- function(object) {
  
  object <- object[object$term == attr(object, "treat"), ]
  object$method_label <- label_ps_method(object$method)
  
  # Plot
  y_lab <- if (attr(object, "exponentiate")) "Hazard ratio" else "Log hazard ratio"
  grp_name <- dplyr::group_vars(object)
  ggplot(object, aes(x = method_label, y = estimate)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), size = .3) +
    facet_wrap(~.data[[grp_name]]) +
    xlab("Method") +
    ylab(y_lab) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1)) 
}