#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

# Generic plotting methods -----------------------------------------------------
#' Plot density
#' 
#' Generic methods for plotting kernel density estimates.
#' 
#' @seealso [plot_density.data.frame()], [plot_density.grouped_psweight_mi()] 
#' @export
plot_density <- function(object, ...) {
  UseMethod("plot_density", object)
}

# Getters ----------------------------------------------------------------------
#' @export
get_ps_catvars <- function() {
   c("race_grouped", "sex", "smoker_grouped", "histology", "stage_grouped")
}

#' @export
get_ps_labs <- function() {
  c("Age", "Days since Dx", "Race", "Sex", "Smoker", "Histology", "Stage")
}

#' @export
get_ps_vars <- function() {
  v <- c("age", "days_since_dx", get_ps_catvars())
  names(v) <- get_ps_labs()
  return(v)
}

#' @export
get_ps_var_list <- function(x) {
  vl <- vector(mode = "list", length = nrow(x)) # List of variables for each analysis
  vi <- 3:(ncol(x) - 1) # Indices for variables 
  for (i in 1:length(vl)) {
    is_v <- x[i, vi] # Vector of booleans indicating whether variable is in model
    v <- names(is_v[which(is_v != 0)]) # Vector of variables actually in model
    n_v <- length(v) # Number of covariates
    vl[[i]] <- rep(NA, n_v)
    for (j in 1:n_v) {
      vl[[i]][j] <- v[j]
    }
  }
  names(vl) <- x$analysis_id
  return(vl)
}

# Labeling ---------------------------------------------------------------------
#' Labels for propensity score methods
#' 
#' Convert a variable containing propensity score methods (as in the argument
#' `methods` from [psweight()]) to a factor variable with pretty labels for 
#' the propensity score methods.
#' 
#' @param x A variable containing the propensity methods.
#' 
#' @return A [`factor`] variable.
#' @examples
#' label_ps_method(c("iptw_att", "iptw_att_trim", 
#'                   "match_nearest", "match_nearest_caliper",
#'                   "match_genetic", "match_genetic_caliper", 
#'                   "ps_stratification", "unadjusted"))
#' @export
label_ps_method <- function(x) {
  levels = c("iptw_att", "iptw_att_trim", 
             "match_nearest", "match_nearest_caliper",
             "match_genetic", "match_genetic_caliper", 
             "ps_stratification", "unadjusted")
  labels <- c("IPTW-ATT", "IPTW-ATT (trim)", 
              "Nearest neighbor matching", "Nearest neighbor matching (caliper)",
              "Genetic matching", "Genetic matching (caliper)", 
              "PS stratification", "Unadjusted")
  index <- match(unique(x), levels)
  
  factor(
    x, 
    level = levels[index],
    labels = labels[index]
  )
}


# Row bind data frame like objects ---------------------------------------------
#' Row bind data frame like objects
#' 
#' Row bind a list of objects that inherit from [`dplyr::tibble`] to create 
#' a single object grouped by id. Objects in the list must currently be of
#' class [`psweight`], [`psweight_mi`], [`surv_ate`], or [`marginal_survival`].
#' 
#' @param l A list containing the objects of interest, which must all be of the
#'  same class and inherit from `dplyr::tibble`.
#' @param id 	Identifier for for the returned object used to group results.
#'  Passed to the `.id` argument of [dplyr::bind_rows()].
#' @param integer_id If `TRUE`, then an intetger sequence is used for `id`; if
#' `FALSE`, then the names of `l` are used.
#'  
#' @return A single object of class `"grouped_{class}"` 
#' where `{class}` is the class of element in `l`. The returned object also 
#' inherits from `dplyr::tibble` and is grouped by `id`.
#' @export
rbind_list <- function(l, id, integer_id = FALSE) {
  class1 <- class(l[[1]])[1]
  if (!class1 %in% c("psweight", "psweight_mi", "surv_ate", "marginal_survival")) {
    stop("Object classes in 'l' are not supported.")
  }
  
  if(any(!sapply(l, inherits, what = class1))) {
    stop(paste0("Each element of 'l' must be of class '", class1, "'"),
         call. = FALSE)
  }
  if (integer_id) names(l) <- NULL
  
  res <- l %>% 
    bind_rows(.id = id) %>%
    group_by(.data[[id]])
  if (integer_id) res[[id]] <- as.integer(res[[id]])
  class(res) <- c(paste0("grouped_", class1), class(res))
  res <- bind_attributes(res, lapply(l, attributes))
  res
}

bind_attributes <- function(x, attr) {
  UseMethod("bind_attributes")
}

bind_attributes.default <- function(x, attr) {
  x
}

bind_attributes.grouped_psweight <- function(x, attr) {
  x_vars <- unique(unlist(lapply(attr, function (z) z[["x_vars"]])))
  attr(x, "x_vars") <- x_vars
  x
}
  
bind_attributes.grouped_psweight_mi <- function(x, attr) {
  bind_attributes.grouped_psweight(x, attr)
}

# Generic HTML table -----------------------------------------------------------
#' HTML table
#' 
#' Create an HTML table from an `R` object using [knitr::kable()] and style it
#' using [kableExtra::kable_styling()].
#' @seealso [html_table.analysis()]
#' @export
html_table <- function (x, ...) {
  UseMethod("html_table", x)
}

#' @rdname html_table
#' @export
html_table.default <- function(x) {
  knitr::kable(x) %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "responsive"))
}