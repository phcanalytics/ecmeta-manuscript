#' Specification of propensity score models.
#'
#' A lookup table used to determine the specification of the propensity score
#' model for each analysis. 
#'
#' @format A `data.frame`. The first two columns are identifiers for an
#'  analysis. The remaining columns are potential
#'  covariates for the propensity score model: a variable is included if it
#'  has a value of 1 and is excluded if it has a value of 0.
#' 
#' @examples 
#' ecmeta.nsclc::psmodel_specification
#' @seealso [add_ps_formula()], [add_impute_formula()]
"psmodel_specification"