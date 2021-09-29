# Helper function --------------------------------------------------------------
#' Interpolate survival curves
#' 
#' Interpolate survival curves to obtain survival probabilities at specific times, 
#' which may differ from the times at which the original curves were evaluated. 
#' The step-function nature of the
#' predictions are preserved so that if `times` lies within `[t1, t2)`, the 
#' survival curves evaluated at `t1` will be used. An efficient approach using
#' `data.table` is employed.
#' @param x Survival curves. Assumed that there are columns for `strata`,
#' `time`, and `surv`.
#' @param times Times at which to subset survival curves to.
#' @param extend Whether to extend survival predictions beyond maximum followup
#' time.
interpolate_surv <- function(x, times, extend = FALSE) {
  # Cut the original times into intervals
  xl <- split(data.table::as.data.table(x), f = x$strata)
  xl <- lapply(xl, function (z) {
    if (!0 %in% z$time) {
      z <- rbind(data.table(strata = z$strata[1], time = 0, surv = 1),
                 z)
    }
    index <- findInterval(times, z$time)
    data.table(strata = z$strata[1],
               time = times,
               surv = z$surv[index])
  })
  y <- data.table::rbindlist(xl)
  
  # Remove "times" if they are beyond what is observed in the survival curves
  # from "x"
  if (extend == FALSE) {
    y <- y[time <= max(x$time)]
  }
  
  # Return
  dplyr::as_tibble(y)
}

# Main method ------------------------------------------------------------------
#' Marginal survival curves
#' 
#' This is a generic method for estimating marginal survival functions by levels
#' of a stratification variable (typically treatment vs. control). 
#' 
#' @param object An object of the appropriate class.
#' @param formula A [`formula`] object as in [`survival::survfit`]. The right
#' hand side of the formula should have a single `strata` term.
#' @param times Unique times to compute survival probabilities at.
#' @param ... Currently unused. 
#' 
#' @return An object of class `marginal_survival`, which is a
#' [`dplyr::tibble`] containing the following columns:
#' 
#' \describe{
#' \item{\{`strata`\}}{The name of the stratification variable. This will be the same
#' as the name of the variable supplied by the user.}
#' \item{`time`}{The time at which survival was evaluated.}
#' \item{`surv`}{The survival probability.}
#' }
#' @export
marginal_survival <- function(object, ...) {
  UseMethod("marginal_survival")
} 

#' @rdname marginal_survival
marginal_survival.psweight <- function(object, formula, times = NULL, ...) {
  environment(formula) <- environment() 
  km <- survival::survfit(formula, data = object, se.fit = FALSE,
                          weights = object$weight)
  marginal_survival(km, times = times)
}

#' @rdname marginal_survival
#' @export
marginal_survival.psweight_mi <- function(object, formula, times = NULL, ...) {
  data_list <- split(object, object$imp)
  data_list <- lapply(data_list, function(x) {
    class(x)[1] <- "psweight"
    x
  })
  
  # Get unique survival times
  if (is.null(times)) {
    surv_object <- eval(parse(text = as.character(formula)[2]), envir = object)
    times <- sort(unique(surv_object[, "time"]))
  }

  # Compute survival probabilities
  surv <- lapply(data_list, function(z) marginal_survival(z, formula = formula, 
                                                          times = times)) %>%
    rbindlist(idcol = "imp")
  surv <- surv[, .(surv = mean(surv)), by = c("strata", "time")]
  surv <- dplyr::as_tibble(surv)
  class(surv) <- c("marginal_survival", class(surv))
  surv
}

#' @rdname marginal_survival
#' @export
marginal_survival.survfit <- function(object, times = NULL, ...) {
  strata_names <- as.integer(gsub(".*=","",names(object$strata)))
  strata_values <- rep(strata_names, times = object$strata)
  
  surv <- dplyr::tibble(
    strata = strata_values,
    time = object$time,
    surv = object$surv
  )
  
  if (!is.null(times)) {
    surv <- interpolate_surv(surv, times = times)
  }
  class(surv) <- c("marginal_survival", class(surv))
  surv
}

# Plotting ---------------------------------------------------------------------
#' Plot of marginal survival curve by group
#' 
#' Quickly plot a marginal survival using [`ggplot2`] with grouped data. Survival
#' curves for each group are stratified by a stratification variable and plots
#' are faceted by group.
#' 
#' @param object An object of class `grouped_marginal_survival`.
#' @param xlab Label for the x-axis.  
#' @param strata_labels Labels for each level of the stratification variable.
#' @param ... Currently unused.
#' 
#' @return A [`ggplot2`] object.
#' 
#' @export
autoplot.grouped_marginal_survival <- function(object, xlab = "Time", 
                                               strata_labels = c("Control", "Treated"),
                                               ...) {
  
  grp_name <- dplyr::group_vars(object)
  
  # Make plot
  ggplot(object, aes(x = time, y = surv, col = factor(strata))) +
    geom_line() +
    facet_wrap(~factor(.data[[grp_name]]), scales = "free_x") +
    xlab(xlab) +
    ylab("Proportion surviving") +
    scale_colour_discrete(name = "", labels = strata_labels) +
    theme(legend.position = "bottom")
  
}