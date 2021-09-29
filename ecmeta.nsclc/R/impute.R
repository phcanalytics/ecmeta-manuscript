# Impute -----------------------------------------------------------------------
#' Formulas for multiple imputation
#' 
#' Create formulas for multiple imputation. `impute_formula()` creates the 
#' imputation formulas as functions of an existing formula object (e.g,
#' a formula used for a propensity score model) and `add_impute_formula()` adds
#' them to an [`analysis`][new_analysis()]
#' object. 
#' @param x An [`analysis`][new_analysis()] object.
#' @param formula A [`formula`] object in `impute_formula` and a column 
#' identifier for the formula objects in `add_impute_formula()`.
#' @param interaction Logical. If `TRUE`, then treatment status is interacted with
#' each covariate, which results in separate imputations models for the trial
#' and external control.
#' @return `impute_formula()` returns a single [`formula`] object.
#' `add_impute_formuls` returns an [`analysis`][new_analysis()] object with a 
#' new list column named `impute_formula`, where each element is the formula used 
#' for imputation of covariates required for a particular pairwise analysis. 
#' 
#' @seealso [`psmodel_specification`]
#' @export
#' 
impute_formula <- function(formula, interaction = TRUE) {
  all_vars <- all.vars(formula)
  
  if (interaction) {
    treat <- all_vars[1]
    all_rhs <- all_vars[-1]
    v <- c(treat, all_rhs, paste0(treat, ":", all_rhs))
    f <- stats::as.formula(paste0("~", paste(v, collapse = "+")))
  } else{
    f <- stats::as.formula(paste0("~", paste(all_vars, collapse = "+")))
  }
  
  return(f)
}

#' @rdname impute_formula
#' @export
add_impute_formula <- function(x, formula, interaction = TRUE) {
  f <- lapply(formula, impute_formula, interaction = interaction)
  x$impute_formula <- f
  x
}

#' Multiple imputation
#' 
#' Multiple imputation of covariates for the propensity score model using 
#' [Hmisc::impute.transcan()]. The model used for multiple imputation must
#' be fit with [Hmisc::aregImpute()].
#' @param object An [`Hmisc::aregImpute`] object.
#' @param data The dataset that will be multiply imputed.
#' @param id Name of the variable in `data` uniquely identifying each row.
#' This is added as a column to the output.
#' @return An "stacked" [`dplyr::tibble`] containing all of the multiply
#' imputed datasets row binded together. A column `imp` identifies the 
#' imputation number.
#' 
#' @export
multi_impute <- function(object, data, id = "patient_id", ...) {
  imputed_data <- vector(mode = "list", length = object$n.impute)
  for (i in 1:length(imputed_data)) {
    imputed_data[[i]] <- as.data.frame(
      Hmisc::impute.transcan(x = object, data = data, 
                             imputation = i, list.out = TRUE, pr = FALSE,
                             ...)
      )
    imputed_data[[i]][[id]] <- data[[id]]
  }
  res <- data.table::rbindlist(imputed_data, idcol = "imp")
  for(v in colnames(res)) {
    if (inherits(res[[v]], "impute"))  data.table::set(res, j = v, value = as.double(res[[v]]))
  } 
  return(dplyr::as_tibble(res))
}

# Check imputations ------------------------------------------------------------
#' Observed vs imputed data
#' 
#' Compare the observed and imputed values of either all of the continuous
#' or all of the categorical variables used in an imputation. 
#' @param object An [`Hmisc::aregImpute`] object.
#' @param data A dataset that was multiply imputed using [multi_impute()].
#' @param type Whether to examine continuous or categorical data. Both cannot
#' be assessed at the same time since they are different data types. 
#' @param labels Optional labels for the variables contained in the columns of
#' `data`.
#' 
#' @return An `observed_vs_imputed` object that inherits from [`dplyr::tibble`]
#' containing the columns:
#' \describe{
#' \item{value}{The value of a particular variable.}
#' \item{type}{Either `"Observed"` or the imputation number (e.g., `"Imputation 1"`,
#' `"Imputation 2"`, ...).}
#' \item{var}{The name of the variable that was imputed.}
#' \item{label}{A pretty label for the variable that was imputed.}
#' }
#' 
#' @export
observed_vs_imputed <- function(object, data, type = c("categorical", "continuous"),
                                labels = NULL) {
  type <- match.arg(type)
  
  # Get variable names
  imputed_varnames <- names(object$nna[object$nna > 0])
  all_varnames <- names(object$imputed)
  cat_levels <- object$cat.levels
  cat_levels[sapply(cat_levels, is.null)] <- NULL
  cat_varnames <- names(cat_levels)
  cont_varnames <- all_varnames[!all_varnames %in% cat_varnames]
  if (type == "continuous") {
    imputed_varnames <- imputed_varnames[imputed_varnames %in% cont_varnames]
  } else {
    imputed_varnames <- imputed_varnames[imputed_varnames %in% cat_varnames]
  }
  
  # Labels for variables
  if (is.null(labels)) {
    imputed_labels <- imputed_varnames
  } else {
    imputed_labels <- names(labels)[match(imputed_varnames, labels)]
  }
  
  # In case where there are no imputed values, return NULL
  if (length(imputed_varnames) == 0) return (NULL)
  
  # Get observed and imputed values for each variable with nonzero missing values
  n_rows <- nrow(data[data$imp == 1, ])
  all_rows <- 1:n_rows
  n_imputed_vars <- length(imputed_varnames)
  res <- vector(mode = "list", length = n_imputed_vars)
  
  for (i in 1:n_imputed_vars) {
    # Get value of variable
    imputed_varname <- imputed_varnames[[i]]
    values <- data[[imputed_varname]]
    
    # Identify all rows with NA in the stacked dataset
    na_rows1 <- object$na[[imputed_varname]] 
    na_rows <- c(sapply(n_rows * 0:(object$n.impute - 1), function (row_start) {
      row_start + na_rows1
    }))
    
    # Get rows with NAs
    imputed_values <- values[na_rows]
    
    # Get observed data
    obs_rows <- all_rows[!all_rows %in% na_rows1]
    observed_values <- data[data$imp == 1, ][[imputed_varname]][obs_rows]
    
    # Combined imputed and observed data
    n_observed <- length(observed_values)
    n_imputed <- length(na_rows1)
    if (type == "categorical") {
      res_value <- c(as.character(observed_values), as.character(imputed_values))
    } else{
      res_value <- c(observed_values, imputed_values)
    }
    
    
    res[[i]] <- dplyr::tibble(
      value = res_value,
      type = rep(c("Observed", paste0("Imputation ", 1:object$n.impute)), 
                   times = c(n_observed, rep(n_imputed, object$n.impute))),
      var = imputed_varname,
      label = imputed_labels[[i]]
    ) 
  }
  res <- dplyr::bind_rows(res)
  class(res) <- c("observed_vs_imputed", class(res))
  return(res)
}

#' Observed vs. imputed plot
#' 
#' A plot comparing observed and imputed data for each variable. There is currently
#' assumed to be a group column (`grp`) which identifies the separate analyses.
#' The plot uses [`ggplot2::facet_grid()`] with groups as rows and variables
#' as columns. Categorical data is plotted with a bar plot and kernel density 
#' estimates are used for the continuous variables.
#' 
#' @param object An [`observed_vs_imputed`] object.
#' 
#' @return A [ggplot2::ggplot] object.
#' @export
autoplot.observed_vs_imputed <- function(object) {
  if (is.factor(object$value) | is.character(object$value)) {
    type <- "categorical"
  } else {
    type <- "continuous"
  }
  
  pdata <- object %>%
    dplyr::mutate(grp = factor(grp, levels = unique(object$grp)))
  
  if (type == "categorical") {
    
    pdata <- pdata %>%
      dplyr::group_by(grp, label, type, value) %>%
      dplyr::summarise(n = dplyr::n()) %>%
      dplyr::group_by(grp, label, type) %>%
      dplyr::mutate(freq = n / sum(n)) 
    
    ggplot(pdata,
           aes(x = value, y = freq, fill = type)) +
      geom_bar(position = "dodge", stat = "identity") +
      facet_grid(grp~label, scales = "free_x", switch = "y") +
      scale_fill_discrete(name = "") +
      xlab("") +
      ylab("Proportion") +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  } else{ # Continuous case
    
    ggplot(pdata,
           aes(x = value, col = type)) +
      geom_density(position = "jitter") +
      facet_grid(grp~label, scales = "free", switch = "y") +
      xlab("") + ylab("Density") +
      scale_color_discrete(name = "") +
      theme(legend.position = "bottom")
  }
}