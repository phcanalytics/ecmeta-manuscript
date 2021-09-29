# Pool estimates across multiple imputations -----------------------------------
pool <- function(q, u) {
  # Estimates
  q <- do.call("rbind", q)
  qbar <- apply(q, 2, mean)
  
  # Covariance matrix
  m <- length(u)
  ubar <- Reduce("+", u) / m # Within imputation variance
  b <- apply(q, 2, var) # Between imputation variance
  t <- ubar + (1 + 1/m) * b
  
  # Return
  return(list(qbar = qbar, t = t))
}

# Estimate Cox models (including within multiple imputation context) -----------
coxph_surv_ate <- function(data, ...) {
  UseMethod("coxph_surv_ate")
} 

#' @importFrom survival strata cluster
coxph_surv_ate.psweight <- function(data, formula, double_adjustment,
                                    n_strata = 10) {
  id <- attr(data, "id")
  x_vars <- attr(data, "x_vars")
  method <- data$method[1]
  
  if (!double_adjustment) {
    fcox <- update(formula, paste("~ . +", "cluster(", id, ")"))
  } else {
    fcox <- update(formula, paste(" ~. +", paste(x_vars, collapse = " + "),
                                  " + cluster(", id, ")"))
  }
  if (method == "ps_stratification") {
    fcox <- update(fcox, paste(" ~. + strata(ps_strata)"))
  }
  
  # Estimation
  environment(fcox) <- environment()
  if (method == "ps_stratification") {
    ps_probs <- seq(0, 1, length.out = n_strata + 1) 
    ps_strata <- cut(data$ps,
                     breaks = quantile(data$ps, probs = ps_probs),
                     include.lowest = TRUE)
  }
  survival::coxph(fcox, data = data, weights = data[["weight"]])
}

coxph_surv_ate.psweight_mi <- function(data, formula, double_adjustment,
                                       n_strata = 10) {
  data_list <- split(data, data$imp)
  lapply(data_list, function (z) {
    class(z)[1] <- "psweight" # Ned to convert from psweight_mi to psweight class
    coxph_surv_ate(z, formula = formula, double_adjustment = double_adjustment,
                   n_strata = n_strata)
  })
}

# Tidy coefficients from Cox models (including within multiple imputation) -----
get_coef <- function(data, ...) {
  UseMethod("get_coef")
} 

get_coef.coxph <- function(object) {
  list(
    estimate = coef(object),
    vcov = vcov(object)
  )
}

tidycoef <- function(data, ...) {
  UseMethod("tidycoef")
} 

tidycoef.default <- function(object) {
  coefs <- get_coef(object)
  
  coef_tbl <- dplyr::tibble(
    term = names(coefs$estimate),
    estimate = coefs$estimate,
    se = sqrt(diag(coefs$vcov)),
  )
  coef_tbl$lower <- coef_tbl$estimate - qnorm(.975) * coef_tbl$se
  coef_tbl$upper <- coef_tbl$estimate + qnorm(.975) * coef_tbl$se
  coef_tbl
}

tidycoef.list <- function(object) {
  # Get coefficients and variance covariance matrix
  n_imputations <- length(object)
  est <- vcov <- vector(mode = "list", length = n_imputations)
  for (i in 1:n_imputations) {
    coefs <- get_coef(object[[i]])
    est[[i]] <- coefs$estimate
    vcov[[i]] <- coefs$vcov
  }
  
  # Pool via Rubin's rule
  coef_pooled <- pool(est, vcov)
  
  # Make tidy table
  coef_tbl <- dplyr::tibble(
    term = names(coef_pooled$qbar),
    estimate = coef_pooled$qbar,
    se = sqrt(diag(coef_pooled$t)),
  )
  coef_tbl$lower <- coef_tbl$estimate - qnorm(.975) * coef_tbl$se
  coef_tbl$upper <- coef_tbl$estimate + qnorm(.975) * coef_tbl$se
  coef_tbl
}

#  Main function to estimate treatment effects ---------------------------------
#' Average treatment effects for survival data
#' 
#' Estimate average treatment effects (ATE) from survival data. A Cox proportional
#' hazards model ([`survival::coxph()`]) is used to estimate hazard ratios and
#' Kaplan-Meier estimators (`survival::survfit.formula()`) are used to estimate 
#' survivor functions. `surv_ate()` estimates ATEs for a single group
#' and `map_surv_ate` iterates over groups and estimates ATEs for each group. 
#' Confidence intervals for hazard ratios are estimated by clustering on `id` using 
#' the `cluster` option in `coxph()`. Estimates and confidences intervals are pooled
#' using Rubin's rule when `object` contains weights from multiply imputed datasets.
#' 
#' @param object,object_list Either an object of class [`psweight`] or [`psweight_mi`] or 
#' a list of such objects.
#' @param ydata,ydata_list Either a data frame containing the response variables 
#' for survival modeling or a list of such objects. Both `object` and `ydata` 
#' must contain the `id` variable.
#' @param response A string containing the response for a survival model.
#' This is the let hand side of a `~` operator.
#' @param id The name of the column in `object` and `ydata` identifying each 
#' patient.
#' @param double_adjustment If `TRUE`, then covariates (in addition to treatment
#' assignment) are used in the Cox models.
#' @param ps_stratification If `TRUE`, then models are fit that stratify on
#' values of the propensity score, with the number of strata determined by 
#' `n_strata`. Default is `FALSE`.
#' @param n_strata Number of strata to divide subjects into based on values of
#' the propensity score.
#' @return A [`dplyr::tibble`] with one row for each `method` in `object`.
#' Each row contains the columns: `method` (the propensity score method),
#' `fit` (`tibble` containing the fitted Cox model),
#' `loghr` (`tibble` containing estimates of the log hazard ratio), and `surv` (A
#' `tibble` containing estimated survival probabilities for both treatment
#' and control).
#' @export
surv_ate <- function(object, ydata, response, id, double_adjustment = FALSE,
                     ps_stratification = FALSE, n_strata = 10) {
  
  # Remove patients with zero weights
  object <- object[object$weight > 0, ]
  
  # Create data for modeling from the initial "object". 
  # There is one dataset for each method. We will create a dummy element
  # for stratification on the propensity score
  data <- dplyr::left_join(object, ydata, by = id) 
  data_list <- split(data, data$method)
  if (ps_stratification) {
    names(data_list)[names(data_list) == "unadjusted"] <- "ps_stratification"
    data_list <- c(data_list, 
                   list(unadjusted = data_list$ps_stratification))
    data_list$ps_stratification$method <- "ps_stratification"
  }
  
  # Estimate treatment effects
  f <- as.formula(paste0(response, "~", attr(object, "treat")))
  
  ## Cox model fits
  fit <- lapply(data_list, function(x) {
    coxph_surv_ate(x, formula = f, double_adjustment = double_adjustment,
                   n_strata = n_strata)
  }) 
  
  loghr <- lapply(fit, function(x) {
    tidycoef(x)
  }) 
  
  ## Survival functions
  if (!double_adjustment) {
    surv <- lapply(data_list, function(x) {
      marginal_survival(x, formula = f)
    })
  }
 
  # Return
  res <- dplyr::tibble(
    method = names(data_list),
    fit = fit,
    loghr = loghr
  )
  if (!double_adjustment) res$surv <- surv
  class(res) <- c("surv_ate", class(res))
  attr(res, "treat") <- attr(object, "treat")
  res
}

#' @rdname surv_ate
#' @export
map_surv_ate <- function(object_list, ydata_list, response, id,
                         double_adjustment = FALSE,
                         grp_id,
                         integer_grp_id = FALSE) {
  args <- match.call(expand.dots = FALSE)
  
  res <- purrr:::map2(object_list, ydata_list, function (x, y){
    surv_ate(x, ydata = y, 
             response = response,
             id = id,
             double_adjustment = double_adjustment)
  }) %>%
    rbind_list(id = grp_id, integer_id = integer_grp_id)
  attr(res, "call") <- match.call()
  res
}

# Survival curves --------------------------------------------------------------
#' @export
plot_survival <- function(object, ...) {
  UseMethod("plot_survival")
}

#' @export
plot_survival.grouped_surv_ate <- function(object, method = "iptw_att_trim",
                                           ...) {
  
  object <- object[object$method == method, ]
  grp_name <- dplyr::group_vars(object)
  
  # Combine survival curves
  names(object$surv) <- object[[grp_name]]
  integer_id <- if (is.integer(object[[grp_name]])) TRUE else FALSE
  surv <- rbind_list(object$surv, id = grp_name, integer_id = integer_id)
  
  # Plot
  autoplot(surv, ...)
}