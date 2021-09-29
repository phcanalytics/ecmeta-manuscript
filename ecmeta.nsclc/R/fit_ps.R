#' Fit propensity score model
#' 
#' Fit a propensity score model using a specified engine. `fit_ps()` fits 
#' a model on a single dataset and `fit_ps_mi()` fits a model on each of 
#' `m` multiply imputed datasets.
#' @param formula An object of class [`formula`] with the treatment variable to the
#' left of the `~` operator and the model terms on the right. 
#' @param engine The "engine" used to fit the propensity score model. Currently
#' limited to logistic regression fit with [stats::glm()].
#' @param data A data frame containing the variables in the model.
#' For `fit_ps_mi()`, there must be an additional column `imp` denoting the
#'  imputation number. 
#' @param ... Additional arguments to pass to [stats::glm()].
#' @return `fit_ps()` returns an object of class `fit_ps` with an element `fit` containing
#' the underlying fitted model. `fit_ps_mi()` returns an object of class`fit_ps_mi`, 
#' which is a list of `fit_ps` objects.
#' 
#' @export
fit_ps <- function(formula, data, engine = "glm", ...) {
  engine <- match.arg(engine)
  if (engine == "glm") {
    fit <- stats::glm(formula, family = binomial(), data = data, ...)
  }
  res <- list(fit = fit)
  class(res) <- "fit_ps"
  return(res)
}

#' @rdname fit_ps
#' @export
fit_ps_mi <- function(formula, data, engine = "glm", ...) {
  if (inherits(data, "data.frame")) {
    data <- split(data, data$imp)
  }
  fits <- lapply(data, function (x) {
    fit_ps(formula, data = x, engine = engine, ...)
  })
  class(fits) <- "fit_ps_mi"
  return(fits)
}