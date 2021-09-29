#' Predict propensity scores
#' 
#' Predict propensity scores from a fitted statistical model. Fitted value from
#' the model fit are used.
#' @param object A fitted propensity score model of the appropriate class. 
#' @param id Name of the variable from the data frame used to fit `object` that 
#' uniquely identifies each row. This is added as a column to the output.
#' @param ... Currently unused. 
#' @return An object of class `ps` or `ps_mi`. A `ps_mi` object is a list
#' of `ps` objects (one for each imputation) and a `ps` object is a list with 
#' elements:
#' \describe{
#' \item{ps}{The propensity score.}
#' \item{lp}{The linear predictor (i.e., the logit of the propensity score).}
#' \item{data}{The data used to fit the propensity score model.}
#' \item{formula}{The model formula used to fit the propensity score model.}
#' \item{treat}{The name of the treatment variable.}
#' \item{x_vars}{The names of all covariates (prior to any transformations or 
#' dummy variable coding) used in the propensity score model.}
#' \item{id}{The }
#' }
#' 
#' @export
predict_ps <- function (object, ...) {
  UseMethod("predict_ps", object)
}

#' @rdname predict_ps
#' @export
predict_ps.fit_ps <- function(object, id = "patient_id", ...) {
  ps <- predict(object$fit, newdata = object$fit$data, type = "response")
  lp <- stats::qlogis(ps)
  all_vars <- all.vars(object$fit$formula)
  
  res <- list(
    ps = ps,
    lp = lp,
    data = object$fit$data,
    formula = object$fit$formula,
    treat = all_vars[1],
    x_vars = all_vars[-1],
    id = id
  )
  class(res) <- "ps"
  return(res)
}

#' @rdname predict_ps
#' @export
predict_ps.fit_ps_mi <- function(object, id = "patient_id", ...) {
  res <- lapply(object, function (x) {
    predict_ps(x, id = id, ...)
  })
  names(res) <- paste0("imputation", 1:length(res))
  class(res) <- "ps_mi"
  return(res)
}