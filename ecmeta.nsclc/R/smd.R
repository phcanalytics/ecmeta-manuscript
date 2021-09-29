# Helper functions -------------------------------------------------------------
reverse_factor <- function(x) {
  factor(x, levels = rev(sort(unique(x))))
}

smd_cov_labels <- function(x) {
  u <- rev(sort(unique(x)))
  u <- c(u[!u %in% c("logit_ps", "ps")], "ps", "logit_ps")
  factor(x, levels = u)
}

# Main function to compute standardized mean differences  ----------------------
#' Standardized mean differences
#' 
#' Compute standardized mean differences (SMDs) between the treated and control
#' subjects. Computed as \eqn{(\mu_t - \mu_c)/\sigma} where \eqn{\mu_t} is the 
#' weighted mean among the treated subjects, \eqn{\mu_c} is the weighted mean
#' among the control subjects, and \eqn{\sigma} is the standard deviation in
#' the unweighted sample (among the treated subjects only since we are interested
#' in the average treatment effect for the treated).
#' @param object An object of the appropriate class.
#' @param x_vars A character vector containing the covariates used to check
#' for balance. These must be columns in `object`. If `NULL`, then variables
#' that were used to fit the propensity score model(s) are used. 
#' @return A `smd` object that inherits from [`dplyr::tibble`]. The following
#' columns are always returned:
#' \describe{
#' \item{var}{A variable - either a covariate or the propensity score.}
#' \item{method}{The propensity score method.}
#' \item{mean_control}{The (weighted) mean among the controls.}
#' \item{mean_treat}{The (weighted) mean among the treated.}
#' \item{sd_treat}{The standard deviation among the treated.}
#' \item{smd}{The standardized mean difference.}
#' }
#' In addition, a column for the group included if `object` is grouped and the
#' column `imp` is included if `object` was instantiated from a propensity
#' score fit using multiply imputed datasets. 
#' @export
smd <- function (object, ...) {
  UseMethod("smd", object)
}

#' @rdname smd
#' @export
smd.grouped_psweight_mi <- function(object, x_vars = NULL) {
  
  grp_name <- dplyr::group_vars(object)
  object <- as.data.table(object)
  
  # Get variables in propensity scores in numeric format
  if (is.null(x_vars)) {
    x_vars <- unique(unlist(attr(object, "x_vars"), use.names = FALSE))
  } 
  f <- stats::as.formula(paste0("~-1+", paste(x_vars, collapse = "+")))
  x <- stats::model.matrix(f, data = model.frame(object, na.action = "na.pass"))
  
  # Reshape into longer format
  id_vars <- c(grp_name, "treat", "imp", "method", "weight")
  z <- cbind(as.data.table(x), object[, c(id_vars, "logit_ps", "ps"), with = FALSE])
  z <- data.table::melt(data.table::as.data.table(z), 
                        id.vars = id_vars,
                        measure.vars = c(colnames(x), "logit_ps", "ps"),
                        variable.name = "var",
                        variable.factor = FALSE,
                        value.name = "value")
  
  # Compute weighted means
  z_wm <- z[, lapply(.SD, weighted.mean, w = weight), 
            by = c("treat", grp_name, "imp", "method", "var"), .SDcols = "value"]
  z_wm[, treat := ifelse(treat == 0, "mean_control", "mean_treat")]
  z_wm <- data.table::dcast(
    z_wm, 
    as.formula(paste0(grp_name, " + imp + method + var ~ treat")),
    value.var = "value"
  )
  
  # Compute standard deviation (in unadjusted treated group for ATT)
  z_sd <-  z[treat == 1 & method == "unadjusted", lapply(.SD, sd), 
             by = c(grp_name, "imp", "var"), .SDcols = "value"]
  data.table::setnames(z_sd, "value", "sd_treat")
  
  # Create final object
  res <- merge(z_wm, z_sd, by = c(grp_name, "imp", "var"))
  res[, smd := (mean_treat - mean_control)/sd_treat]
  
  # Return
  res <- dplyr::as_tibble(res) %>%
    dplyr::group_by(.data[[grp_name]])
  class(res) <- c("smd", class(res))
  return(res)
}

#' Standardized mean difference list
#' 
#' A list of objects of class [`smd`].
#' @param ... Objects, possibly named.
#' @return An object of class `smd_list`.
#' @export
smd_list <- function(...) {
  objects <- list(...)
  class(objects) <- "smd_list"
  return(objects)
}

# Plotting functions -----------------------------------------------------------
#' Plot standardized mean differences
#' 
#' Quickly plot standardized mean differences. Plots are faceted by group. 
#' @inheritParams plot_ps
#' @param id Character vector of the same length as `object` to name each element
#' in the list of `smd` objects. Passed to the `.id` argument in 
#' `dplyr::bind_rows()`.
#' @details If `method` is `NULL`, then all methods are compared in a single plot.
#' However, in this case it is infeasible to display each covariate in the plot
#' so only the (logit of the) propensity scores are compared. If `object` is of
#' class `smd_list`, then SMDs are compared (via the `colour` aesthetic) across 
#' all objects in the list.
#' @return A [`ggplot2::ggplot`] object.
#' @export
autoplot.smd <- function(object, method = "iptw_att") {
  grp_name <- dplyr::group_vars(object)
  
  # Subset according to method argument
  if (!is.null(method)) {
    pdata <- object[object$method == method, ]
  } else {
    pdata <- object[object$var == "logit_ps", ]
  }
  
  # Summarize data for plotting
  pdata <- pdata %>%
    dplyr::filter(!is.na(smd)) %>%
    dplyr::mutate(abs_smd = abs(smd)) %>%
    dplyr::group_by(.data[[grp_name]], var, method) %>%
    dplyr::summarise(mean_abs_smd = mean(abs_smd),
                     min_abs_smd = min(abs_smd),
                     max_abs_smd = max(abs_smd)) 
  
  # Plotting
  add_layers <- function(p) {
    p + geom_point(position = position_dodge(width = 0.3)) +
      geom_linerange(aes(xmin = min_abs_smd, xmax = max_abs_smd),
                     position = position_dodge(width = 0.3)) +
      geom_vline(xintercept = .1, linetype = "dashed", color = "grey") +
      facet_wrap(~.data[[grp_name]], scales = "fixed") +
      ylab("") +
      theme(legend.position = "none")    
  }
  
  if (!is.null(method)) { # Plot for single method
    pdata$var <- smd_cov_labels(pdata$var)
    p <- ggplot(pdata,
                aes(x = mean_abs_smd, y = var)) %>%
      add_layers() +
      xlab("Absolute standardized mean difference") 
  } else { # Plot summarizing across methods
    pdata$method <- reverse_factor(label_ps_method(pdata$method))
    p <- ggplot(pdata,
                aes(x = mean_abs_smd, y = method, col = method)) %>%
      add_layers() +
      xlab("Absolute standardized mean difference of the logit of the PS") 
  }
  
  return(p)
}

#' @rdname autoplot.smd
#' @export
autoplot.smd_list <- function(object, method = "iptw_att", id = NULL) {
  
  grp_name <- dplyr::group_vars(object[[1]])
  
  # Combine lists
  class(object) <- "list"
  if (!is.null(id)) names(object) <- id
  pdata <- dplyr::bind_rows(object, .id = "id")
  
  # Subset to method of choice
  pdata <- pdata[pdata$method == method, ]
  
  # Summarize data for plotting
  pdata <- pdata %>%
    dplyr::filter(!is.na(smd)) %>%
    dplyr::mutate(abs_smd = abs(smd)) %>%
    dplyr::group_by(id, .data[[grp_name]], var, method) %>%
    dplyr::summarise(mean_abs_smd = mean(abs_smd),
                     min_abs_smd = min(abs_smd),
                     max_abs_smd = max(abs_smd)) 
  
  # Plot
  pdata$var <- smd_cov_labels(pdata$var)
  p <- ggplot(pdata,
              aes(x = mean_abs_smd, y = var, col = id)) +
    geom_point(position = position_dodge(width = 0.3)) +
    geom_linerange(aes(xmin = min_abs_smd, xmax = max_abs_smd),
                   position = position_dodge(width = 0.3)) +
    geom_vline(xintercept = .1, linetype = "dashed", color = "grey") +
    facet_wrap(~.data[[grp_name]], scales = "fixed") +
    xlab("Absolute standardized mean difference") +
    ylab("") +
    scale_color_discrete(name = "") +
    theme(legend.position = "bottom")    
  
  return(p)
  
}

