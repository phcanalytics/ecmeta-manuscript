# Main methods to compute weights  ---------------------------------------------
#' Propensity score weighting
#' 
#' A generic function used to weight data using estimated propensity scores. Six
#' separate propensity score methods can be employed: (i) inverse probability weighting
#' for the average treatment effect on the treated (IPTW-ATT weighting), (ii)
#' IPTW-ATT weighting where controls with predicted propensity scores less than
#' the 1st percentile or greater than the 99th percentile are trimmed (i.e., excluded), 
#' (iii) 1:1 nearest neighbor matching without a caliper, (iv) 1:1 nearest 
#' neighbor matching with a caliper (of 0.25 standard deviations of the
#' linear propensity score), (v) 1:1 genetic matching without a caliper,
#' and (vi) 1:1 genetic matching with a caliper. Matching is always performed
#' without replacement unless the number of controls is less then the number of 
#' treated (in which case it is performed with replacement). Unadjusted 
#' (i.e., unweighted) results are also always returned.
#' 
#' @param object An object of the appropriate class containing the propensity
#' score to match on and other relevant information.
#' @param methods The propensity score methods to use. Can select any combination
#' of the available methods. The default is to use IPTW-ATT weighting (with
#' and without trimming) and nearest neighbor matching (with and without a 
#' caliper). Genetic matching can be employed using the options `match_genetic`
#' and `match_genetic_caliper`.
#' @param progress Logical. If `TRUE`, the the iteration number from the list
#' is printed to the screen and save to the outfile.
#' @param outfile File path when `progress = TRUE`.
#' @param ... Additional arguments to pass to [Matching::GenMatch()].
#' 
#' @return A `dplyr::tibble` is always returned that contains (i) each
#' variable included in `data`, (ii) the columns `ps` and `logit_ps` and for the 
#' propensity score and logit of the propensity score, respectively,
#' (iii) the column `weight` for the weight generated, (iv) the column `method`
#' for the propensity score method used, (v) the column `exclude` denoting
#' whether patients were excluded using either trimming or a caliper, and (vi)
#' the column `method_type` for the type of propensity score method (weighting,
#' matching, or none).
#' 
#' The returned object always inherits from [`dplyr::tibble`] but the class
#' depends on `object`. If `object` is of class [`ps_fit`], then there is one row 
#' for each patient and the returned object is of class [`psweight`].  If 
#' `object` is of class [`ps_fit_mi`], then an additional column `imp` is added
#' for the imputation number and an object of class `psweight_mi` is returned. Finally,
#' if `object` is of class [`list`], then an additional column `grp` is added and
#' an object of class `grouped_psweight` or `grouped_psweight_mi` is returned. 
#' 
#' @export
psweight <- function (object, ...) {
  UseMethod("psweight", object)
}

#' @rdname psweight
#' @export
psweight.ps <- function(object, methods = c("iptw_att", "iptw_att_trim",
                                            "match_nearest", "match_nearest_caliper"), ...) {
  method_choices <- c("iptw_att", "iptw_att_trim", 
                      "match_nearest", "match_nearest_caliper", 
                      "match_genetic","match_genetic_caliper")
  methods <- match.arg(methods, method_choices,  several.ok = TRUE)
  methods <- c(methods, "unadjusted")
  
  # Generate x matrix if using genetic matching
  if (any(methods %in% c("match_genetic", "match_genetic_caliper"))) {
    x <- cbind( # Variables to match on in addition to PS
      lp = object$lp,
      stats::model.matrix(object$formula, data = object$data)[, -1]
    )
  } else {
    x <- NULL
  }
  
  # Match with or without replacement? Only replace if number of controls is 
  # less than the number of treated
  if (any(methods %in% c("match_nearest", "match_nearest_caliper",
                         "match_genetic", "match_genetic_caliper"))) {
    n_treat <- table(object$data[[object$treat]])
    replace <- if (n_treat["1"] > n_treat["0"]) TRUE else FALSE
  } else {
    replace <- NULL
  }
  
  # Get weights
  w_list <- lapply(methods, get_weight, object = object, x = x, replace = replace,
                   ...)
  
  # Return
  n_methods <- length(methods)
  
  rep_df <- function(x, ...) {
    x[rep(1:nrow(x), ...), ]
  }
  res <- dplyr::tibble(
    rep_df(object$data, times = n_methods),
    logit_ps = rep(object$lp, times = n_methods),
    ps = rep(object$ps, times = n_methods),
    weight = unlist(w_list),
    method = rep(methods, each = nrow(object$data))
  )
  
  class(res) <- c("psweight", class(res))
  attr(res, "treat") <- object$treat
  attr(res, "x_vars") <- object$x_vars
  attr(res, "id") <- object$id
  return(res)
}

#' @rdname psweight
#' @export
psweight.ps_mi <- function(object, ...) {
  res <- lapply(object, function (x) {
    psweight(x, ...)
  })
  res <- dplyr::bind_rows(res)
  class(res)[1] <- "psweight_mi"
  return(res)
}

#' @rdname psweight
#' @export
psweight.list <- function(object, progress = FALSE, outfile = "psweight_out",
                          ...) {
  i <- 1
  if (progress) cat("Starting propensity score weighting:\n", 
                    file = outfile)
  res <- lapply(object, function (x) {
    out <- psweight(x, ...)
    if (progress) cat(paste0("Completed list element ", i, "\n"), 
                      file = outfile, append = TRUE)
    i <<- i + 1
    return(out)
  })
  return(res)
}

# Functions to get weights for each PS method ----------------------------------
get_weight <- function (object, method, x = NULL, replace = NULL, ...) {
  # Inputs required to generate weights
  ps <- object$ps
  lp <- object$lp
  treat <- object$data[[object$treat]]
  
  # Weight based on method
  if (method == "iptw_att") {
    w <- get_weight_iptw_att(ps, treat = treat)
  } else if (method == "iptw_att_trim") {
    w <- get_weight_iptw_att(ps, treat = treat,
                             trim = TRUE)
  } else if (method == "match_nearest") {
    w <- get_weight_match_nearest(lp, treat = treat, caliper = NULL,
                                  replace = replace)
  } else if (method == "match_nearest_caliper") {
    w <- get_weight_match_nearest(lp, treat = treat, caliper = .25,
                                  replace = replace)
  } else if (method == "match_genetic") {
    w <- get_weight_match_genetic(treat = treat, X = x, 
                                  caliper = NULL, replace = replace, ...)
  } else if (method == "match_genetic_caliper") {
    w <- get_weight_match_genetic(treat = treat, X = x, 
                                  caliper = .25, replace = replace, ...)
  } else if (method == "unadjusted") {
    w <- rep(1, length(treat))
  }
  return(w)
}

trim_iptw_att <- function(w, ps, treat, prob = .99) {
  ps_upper <- quantile(ps, prob)
  ps_lower <- quantile(ps, 1 - prob)
  w <- ifelse(treat == 0 & (ps > ps_upper | ps < ps_lower), 0, w) # Only exclude controls
  return(w)
}

get_weight_iptw_att <- function(ps, treat, trim = FALSE) {
  w <- ifelse(treat == 1, 1,  ps/(1 - ps))
  if (trim) w <- trim_iptw_att(w, ps, treat)
  return(w)
}

get_weight_Match <- function(x) {
  # Get weights
  w <- rep(0, x$orig.nobs)
  
  ## All treated subjects receive a weight of 1
  w[x$index.treated] <- 1 
  
  ## All dropped subjects receive a weight of 0
  w[x$index.dropped] <- 0 
  
  ## Sum weights for each control subject
  w_control <- dplyr::tibble(index = x$index.control, weight = x$weights) %>%
    dplyr::group_by(index) %>%
    dplyr::summarise(weight = sum(weight))
  w[w_control$index] <- w_control$weight
  
  # Return
  return(w)
}

get_weight_match_nearest <- function(lp, treat, M = 1, caliper = .25, 
                                     replace = FALSE) {
  
  m <- Matching::Match(
    Tr = treat,
    X = lp, # Linear propensity score
    M = M,
    caliper = caliper,
    ties = FALSE,
    replace = replace,
    version = "fast",
    estimand = "ATT"
  )
  
  return(get_weight_Match(m))
}

get_weight_match_genetic <- function(treat, X, M = 1, caliper = .25, 
                                     replace = FALSE, ...) {
  
  # Get weights for each covariate
  g <- Matching::GenMatch(
    Tr = treat,
    X = X, # Linear propensity score + other variables to match on
    M = M,
    caliper = caliper,
    ties = FALSE,
    replace = replace,
    estimand = "ATT",
    ...
  )
  
  # Use weights in nearest neighbor matching
  m <- Matching::Match(
    Tr = treat,
    X = X, # Same as above
    M = M,
    caliper = caliper,
    ties = FALSE,
    replace = replace,
    version = "fast",
    estimand = "ATT",
    Weight.matrix = g
  )
  
  return(get_weight_Match(m))
}

# Plot distribution of estimated propensity score weights ----------------------
#' Plot propensity score weights
#' 
#' A histogram displaying propensity score weights. If `object` is of class 
#' `grouped_psweight`, then histograms are faceted by the group number. If
#' `object` is of class `grouped_psweight_mi`, then weights from the imputations
#' are pooled. 
#' @param object An object of the appropriate class.
#' @param method A character vector of length one indicating the propensity score
#'  method used to generate the weights.
#' @return A [`ggplot2::ggplot`] object.

#' @export
plot_weights <- function (object, ...) {
  UseMethod("plot_weights", object)
}

#' @rdname plot_weights
#' @export
plot_weights.grouped_psweight_mi <- function(object, method = "iptw_att") {
  if (length(method) > 1) stop("'method' must be of length 1.")
  
  # Subset data
  ## Don't need unadjusted weights
  pdata <- object[object$method != "unadjusted", ]
  
  ## External controls only since using ATT weights
  treat <- attr(object, "treat")
  pdata <- pdata[pdata[[treat]] == 0, ]
  
  ## By method
  pdata <- pdata[pdata$method == method, ]
  
  # Plot histograms
  grp_name <- dplyr::group_vars(object)
  ggplot(data = pdata, aes(x = weight)) +
    geom_histogram(alpha = .8, bins = 30, colour = "white") +
    xlab("Weights") +
    ylab("Count")  +
    facet_wrap(~pdata[[grp_name]], scales = "free") +
    theme(legend.position = "bottom")
}

# Plot distribution of the propensity score ------------------------------------
#' Plot propensity scores
#' 
#' Kernel density estimates of the propensity scores by treatment assignment in 
#' order to check whether the estimated propensity scores are well balanced. 
#' If `object` is of class `grouped_psweight`, then density plots are faceted by 
#' the group number. If `object` is of class `grouped_psweight_mi`, then propensity
#' scores from the imputations are pooled. 
#' @inheritParams plot_weights
#' @param logit Logical. If `TRUE`, then the logit of the propensity score is 
#' plotted.
#' @return A [`ggplot2::ggplot`] object.
#' @export
plot_ps <- function (object, ...) {
  UseMethod("plot_ps", object)
}

#' @rdname plot_ps
#' @export
plot_ps.grouped_psweight_mi <- function(object, method = "iptw_att", 
                                        logit = FALSE) {
  # Subset data
  treat <- attr(object, "treat")
  pdata <- object[object$method == method, ]
  
  # Plot histograms
  if (!logit) {
    ps_var <- "ps"
    ps_lab <- "Propensity score"
  } else {
    ps_var <- "logit_ps"
    ps_lab <- "Logit of propensity score"
  }
  
  grp_name <- dplyr::group_vars(object)
  ps_var <- if (logit) "logit_ps" else "ps"
  ggplot(data = pdata, aes(x = .data[[ps_var]], col = factor(treat))) +
    geom_density(aes(weight = weight)) +
    scale_color_discrete(name = "", labels = c("Control", "Treated")) +
    xlab(ps_lab) +
    ylab("Density")  +
    facet_wrap(~pdata[[grp_name]], scales = "free") +
    theme(legend.position = "bottom") 
}

# Density plots of continuous covariates ---------------------------------------
#' Density plots to check balance
#' 
#' Check the balance of continuous covariates with density plots.
#' @export
plot_density.grouped_psweight_mi <- function(object, method = "iptw_att", 
                                             var, xlab = NULL) {
  # Subset data
  grp_name <- dplyr::group_vars(object)
  treat <- attr(object, "treat")
  method_str <- method
  pdata <- object %>%
    dplyr::filter(method == method_str) %>%
    dplyr::select(dplyr::all_of(c(var, "weight", grp_name, treat)))
  
  # Plot
  if (is.null(xlab)) xlab <- var
  ggplot(data = pdata, aes(x = .data[[var]], col = factor(treat))) +
    geom_density(aes(weight = weight)) +
    scale_color_discrete(name = "", labels = c("Control", "Treated")) +
    xlab(xlab) +
    ylab("Density")  +
    facet_wrap(~.data[[grp_name]], scales = "free") +
    theme(legend.position = "bottom") 
}