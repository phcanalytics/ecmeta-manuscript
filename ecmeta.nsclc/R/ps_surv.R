#' Propensity score survival analysis
#' 
#' Run a complete pipeline required to estimate average treatment effects (ATEs) with
#' survival data using propensity score methods. The pipeline 
#' (1) fits the propensity score model with [fit_ps()] or [fit_ps_mi()], 
#' (2) predict the propensity score with [predict_ps()], (3) weights
#' patients using the propensity scores with [psweight()], and estimates
#' ATEs with [surv_ate()]. `ps_surv()` implements the pipeline for a single group
#' while `map_ps_surv()` iterates over multiple groups. 
#'
#' @param xdata,xdata_list A dataset containing the variables for propensity
#' score modeling. This can either be an object produced by `make_xdata()`
#' or a standard data frame. If a standard data frame, then `formuls` cannot be
#' `NULL` and `make_xdata()` is run using `formula`.
#' @param formula The propensity score formula as described in [fit_ps()]. Only
#' required if `data` was not produced by `make_xdata()`.
#' @param ydata,ydata_list A data frame or list of data frames containing the
#' response variables as described in [surv_ate()].
#' @param  A string containing the response for a survival model as described in
#' [surv_ate()].
#' @param id The name of the column in `xdata` and `ydata` identifying each patient.
#' @param ... Additional arguments to pass to `make_xdata()`.
#' 
#' @return `ps_surv()` returns an object of class `ps_surv` and `map_ps_surv()`
#' returnes an object of class `grouped_ps_surv`. Both are lists containing
#' the elements:
#' \describe{
#' \item{psfit}{The propensity score model fit using [fit_ps()] or [fit_ps_mi()].}
#' \item{psweight}{The propensity score weights returned by [psweight()].}
#' \item{ate}{The average treatment effects returned by [surv_ate()].}
#' \item{call}{The function call returned by [match.call()].}
#' }
#' @export
ps_surv <- function(xdata, formula = NULL, ydata, response, id = "patient_id",
                    ...) {
  
  if (!inherits(xdata, "model_data")) {
    xdata <- make_xdata(formula = formula, data = xdata, ... )
  } 
  ps_formula <- attr(xdata, "formula")
  
  if (attr(xdata, "imputation") == "multi") {
    ps_fit <- fit_ps_mi(ps_formula, data = xdata)
  } else {
    ps_fit <- fit_ps(ps_formula, data = xdata)
  }
  
  ps <- predict_ps(ps_fit, id = id)
  psw <- psweight(ps)
  ate <- surv_ate(psw, ydata, response = response, id = id)
  res <- list(
    psfit = ps_fit,
    psweight = psw,
    ate = ate,
    call = match.call()
  )
  class(res) <- "ps_surv"
  res
}

#' @rdname ps_surv
#' @export
map_ps_surv <- function(xdata_list, formula = NULL, ydata_list, 
                         response, id = "patient_id",
                         grp_id, integer_grp_id = FALSE, ...) {

  
  if (is.null(formula)) { # Empty list of NULLs
    formula <- vector(mode = "list", length = length(xdata_list))
  }
  res <- purrr::pmap(list(xdata_list, ydata_list, formula), function (x, y, z){
    ps_surv(xdata = x, ydata = y, formula = z, response = response,
            id = id)
  }) 
  
  # Combine list items and return
  psfit <- lapply(res, function (z) z[["psfit"]])
  psw <- lapply(res, function (z) z[["psweight"]]) %>%
    rbind_list(id = grp_id, integer_id = integer_grp_id)
  ate <- lapply(res, function (z) z[["ate"]]) %>%
    rbind_list(id = grp_id, integer_id = integer_grp_id)
  
  
  res <- list(
    psit = psfit,
    psw = psw,
    ate = ate,
    call = match.call()
  )
  class(res) <- "grouped_ps_surv"
  res
}