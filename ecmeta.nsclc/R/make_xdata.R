#' Model data
#' 
#' Create data containing the treatment variable and all other covariates. This 
#' data is required to build the design (i.e., `x`) matrix in either a propensity
#' score model or a covariate adjusted outcomes model. Only variables required
#' by the model are kept. 
#' 
#' @param formula A formula for the propensity score with the treatment variable
#' to the left of the `~` operator and the model terms on the right.
#' @param data The dataset containing the variables specified in `formula`.
#' @param imputation Imputation strategy to use. Use `"multi"` to multiply impute
#' data with [multi_impute()]. No imputation is performed when the `"none"` is 
#' chosen.
#' @param missing_dummy Whether to create dummy variables for variables
#' when observations are missing.
#' @param id The variable uniquely identifying each patient.
#' @param n_impute Number of imnputations to perform. Only used when 
#' `imputation = "multi"`.
#' 
#' @return An object of class `model_data` inheriting from [`dplyr::tibble`]. 
#' It contains all variables used in `formula` and, if multiple imputation
#' is used, an additional column named `imp` indicating the imputation
#' number. 
#' @export
make_xdata <- function(formula, data, 
                       imputation = c("multi", "none"),
                       missing_dummy = FALSE,
                       id = "patient_id",
                       n_impute = 5) {
  imputation <- match.arg(imputation)
  
  # Only use variables needed by formula
  model_data <- stats::get_all_vars(formula, data)
  
  # A bit of variable checking
  for (v in colnames(model_data)) {
    if (is.factor(model_data[[v]])) model_data[[v]] <- droplevels(model_data[[v]])
    if (is.character(model_data[[v]])) model_data[[v]] <- as.factor(model_data[[v]])
  }
  model_data <- bind_cols(model_data, data[, id])
  
  # Add missing indicator categories
  na_dummies <- NULL
  if (missing_dummy) {
    has_na <- sapply(model_data, function (v) sum(is.na(v)) > 0)
    
    make_missing_dummies <- function(x, vars) {
      if (length(vars) == 0) return(NULL)
      m <- matrix(NA, nrow = nrow(x), ncol = length(vars))
      colnames(m) <- vars
      for (v in vars) {
        m[, v] <- ifelse(is.na(x[[v]]), 1, 0)
      }
      m
    }
    
    ## Continuous dummies
    is_cont <- sapply(model_data, is.numeric)
    has_na_cont <- has_na & is_cont
    na_cont <- names(has_na_cont)[has_na_cont]
    na_cont_dummies <- make_missing_dummies(model_data, na_cont)
    
    ### If we are not imputing, then we need to add a constant to the missing
    ### values
    if (imputation == "none") {
      for (v in na_cont) model_data[[v]] <- ifelse(is.na(model_data[[v]]), 0, model_data[[v]])
    }
    
    ## Factor dummies
    is_factor <- sapply(model_data, is.factor)
    has_na_factor <- has_na & is_factor
    na_factor <- names(has_na_factor)[has_na_factor]
    
    ### If we are not imputing, then we add a "NA" level to the existing factors
    ### variables. Otherwise we create the dummies manually
    if (imputation == "none") {
      for (v in na_factor) model_data[[v]] <- addNA(model_data[[v]])
    } else {
      na_factor_dummies <- make_missing_dummies(model_data, na_factor)
    }
    
    ## Combine continuous and factor dummies
    if (imputation == "none") {
      na_dummies <- na_cont_dummies
    } else {
      na_dummies <- cbind(na_factor_dummies, na_cont_dummies)
    }
    if (!is.null(na_dummies)) colnames(na_dummies) <- paste0(".missing_", colnames(na_dummies))
  } # End missing_dummy if statement
                                        
  # Imputation
  if (imputation == "multi") {
    impute_form <- impute_formula(formula)
    impute_fit <- Hmisc::aregImpute(impute_form, data = model_data, 
                                    n.impute = n_impute,
                                    nk = 3, match = "closest", 
                                    tlinear = TRUE, pr = FALSE)
    model_data <- multi_impute(impute_fit, data = model_data) 
  }
    
  # Combine model data and dummies
  if (imputation == "multi" & !is.null(na_dummies)) {
    model_data <- cbind(
      model_data,
      na_dummies[rep(1:nrow(na_dummies), times = n_impute), , drop = FALSE]
    )
  } else if (!is.null(na_dummies)) {
    model_data <- cbind(model_data, na_dummies)
  }
  
  # Update model formula to include missing dummies if needed
  if (!is.null(na_dummies)) {
    formula <- update(
      formula, 
      paste("~ . +", paste(colnames(na_dummies), collapse = " +"))
    )
  }
  
  # Return
  model_data <- dplyr::as_tibble(model_data)
  attr(model_data, "formula") <- formula
  attr(model_data, "imputation") <- imputation
  attr(model_data, "id") <- id
  class(model_data) <- c("model_data", class(model_data))
  model_data
}