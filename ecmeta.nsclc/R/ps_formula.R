#' Propensity score formulas
#' 
#' Formulas for modeling the propensity score. `get_ps_formula()` obtains 
#' the formulas and `add_ps_formula()` adds them to an existing 
#' [`analysis`][new_analysis()] object. 
#' @param x An [`analysis`][new_analysis()] object.
#' @param splines Whether to use natural cubic splines ([splines::ns()]) for
#' continuous variables.
#' @param interaction Whether to include an interaction between cancer stage
#' and days since diagnosis. 
#' @param treat Name of the treatment variable. Default is `"treat"`.
#' @param keep_names If `TRUE`, then `analysis_id` is used to name the list. 
#' Otherwise an unnamed list is returned.
#' 
#' @return `get_ps_formula()` returns a list of formulas. `add_ps_formula()`
#'  returns an [`analysis`][new_analysis()] object 
#' with a new list column named `ps_formula` where each element is the formula
#'  used for a particular pairwise analysis. 
#' 
#' @seealso [`psmodel_specification`]
#' @importFrom splines ns
#' @export
#' @name ps_formula
get_ps_formula <- function(splines = TRUE, interaction = TRUE,
                           treat = "treat", keep_names = FALSE) {
  
  # Read in information required to create formula
  psf_tbl <- ecmeta.nsclc::psmodel_specification
  
  # If interaction = FALSE, don't allow interactions
  if (!interaction) psf_tbl$interaction <- 0
  
  # Create list containing variables in each analysis
  vl <- get_ps_var_list(psf_tbl)
  
  # Convert to a list of formulas
  ## Helper function for a single row
  as_formula <- function(z, interaction) {
    # Splines for continuous variables
    cont_vars <- c("age", "days_since_dx")
    which_days_since_dx <- which(z == "days_since_dx")
    for (j in 1:length(z)) {
      if (z[j] %in% cont_vars && splines) z[j] <- paste0("ns(", z[j], ", df = 5)")
    }
    
    # Interaction terms
    if (interaction) {
      z[length(z) + 1] <- paste0("stage_grouped", ":", z[which_days_since_dx])
    }
    
    # Convert to formula and return
    fz <- stats::as.formula(paste0(treat, "~", paste(z, collapse = "+")))
    return(fz)
  }
  f <- mapply(as_formula, vl, psf_tbl$interaction)
  if (!keep_names) names(f) <- NULL
  return(f)
}
 
#' @rdname ps_formula
#' @export
add_ps_formula <- function(x, splines = TRUE, interaction = TRUE) {

  f <- get_ps_formula(splines = splines, interaction = interaction,
                      keep_names = TRUE)
  f <- f[order(match(names(f), x$analysis_id))] # Reorder based on x if needed
  names(f) <- 1:length(f)
  x$ps_formula <- f
  return(x)
}
