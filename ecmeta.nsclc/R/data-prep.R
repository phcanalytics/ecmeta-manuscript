# Initial preprocessing of NSCLC external control datasets ---------------------
#' Preprocess data
#'
#' Preprocess the external control datasets for propensity score modeling. This 
#' involves combining (i.e., grouping) categorical variables, deriving continuous
#' variables (i.e., time since initial diagnosis), and identifying each of the 
#' distinct analyses for the study (defined by a unique clinical trial and 
#' pairwise comparison between an experimental and comparator arm).
#' @param data The unprocessed [`ecdata::nsclc`] dataset.
#' @param combine_levels If `TRUE`, then levels of factor variables with
#' low counts (typically less than 10) are combined into a single category.
#' This current only applies to the race variable.
#' @param drop_nos In some instance, patients with a not otherwise specified 
#' (NOS) histology are part of the external control but not part of the trial.
#' If `TRUE`, then these patients are removed from the external control cohort;
#' if `FALSE`, then they are not.
#' @return A [`dplyr::tibble`].
#' @export
preprocess <- function(data, combine_levels = TRUE, drop_nos = TRUE) {
  
  # Clean histology 
  ## Only use for 2 studies since it was part of the I/E for the others
  data <- data %>% 
    dplyr::mutate(
      histology = ifelse(!dataset_id %in% c("nct02008227_fi", "nct01903993_fi"),
                       NA, histology) 
    ) 
  
  ## There are no NOS patients in the trials (only in the EC arm)
  if (drop_nos) {
    data <- data %>%
      dplyr::filter(!(dataset_id %in% c("nct02008227_fi", "nct01903993_fi") &
                        histology == "NOS")) 
  }
  
  # Group race
  data <- data %>% dplyr::mutate(
    race_grouped = dplyr::case_when(
      race == "White" ~ "White",
      race == "Black or African American" & dataset_id == "nct02367781_fi" ~ "Black",
      race == "Asian" & dataset_id == "nct01351415_fi" ~ "Asian",
      race != "Asian" & dataset_id == "nct01351415_fi" ~ "Non-Asian",
      is.na(race) ~ NA_character_,
      TRUE ~ "Other" # This should never happen
    )
  )
  
  # Group smoking status
  data <-  data %>% dplyr::mutate(
    smoker_grouped = dplyr::case_when(
      smoker %in% c("Former", "former", "Current", "Current/former") ~ "Current/former",
      smoker == "Never" ~ "Never",
      is.na(smoker) | dataset_id == "nct01493843_fi" | smoker == "" ~  NA_character_,
      TRUE ~ smoker
    )
  )
  
  # Group cancer stage
  data <- data %>% dplyr::mutate(
    stage_grouped = dplyr::case_when(
      stage %in% c("0", "I", "IA", "IB", "II", "IIA", "IIB", "III", "IIIA") ~ "Early",
      stage %in% c("IIIB", "IIIC", "IV", "IVA", "IVB") ~ "Advanced",
      is.na(stage) ~ NA_character_,
      TRUE ~ stage
    )
  )
  
  # Combine factor categories with small sample sizes (only race)
  if (combine_levels) {
    data <- data %>% dplyr::mutate(
      race_grouped = dplyr::case_when(
        race_grouped == "Black" & dataset_id %in% c("nct02008227_fi", "nct01903993_fi",
                                                    "nct02366143_fi", "nct01519804_fi",
                                                    "nct01496742_fi", "nct01366131_fi",
                                                    "nct01493843_fi", "nct02657434_fi") ~ "Other",
        TRUE ~ race_grouped
      )
    )
  }
  
  # Time since initial diagnosis
  data <- data %>% dplyr::mutate(
    days_since_dx = dplyr::case_when(
      source_type == "RCT" & dataset_id != "nct01493843_fi" ~ rand_date - dx_date,
      source_type == "RCT" & dataset_id == "nct01493843_fi" ~ rand_date - metastatic_dx_date,
      source_type == "External" ~ tx_start_date - dx_date,
      TRUE ~ NA_real_
    ),
    days_since_dx = as.numeric(days_since_dx),
    days_since_dx_pos = ifelse(days_since_dx <= 0, 1/365.25, days_since_dx),
    log_days_since_dx = log(days_since_dx_pos),
    days_since_dx_pos = NULL
  ) 
  
  # Clean arm names
  data <- data %>% 
    dplyr::mutate(arm = ifelse(arm %in% c("Carboplatin,Pemetrexed", "Cisplatin,Pemetrexed") &
                               dataset_id == "nct02657434_fi",
                               "Carboplatin or Cisplatin + Pemetrexed",
                               arm))
  
  # Treatment indicator
  data <- data %>%
    dplyr::mutate(treat = ifelse(arm_type == "Experimental", 1, 0))
  
  # IDs for pairwise comparisons
  data <- data %>% dplyr::mutate(
    analysis_num = dplyr::case_when(
      dataset_id == "nct02008227_fi" ~ 1L, # OAK
      dataset_id == "nct01903993_fi" ~ 2L, # POPLAR
      dataset_id == "nct02366143_fi" ~ 3L, # IMpower150
      dataset_id == "nct01351415_fi" ~ 4L,
      dataset_id == "nct01519804_fi" ~ 5L,
      dataset_id == "nct01496742_fi" &
        arm %in% c("Placebo + Bevacizumab + Platinum",
                   "MetMAb + Bevacizumab + Platinum + Paclitaxel") ~ 6L,
      dataset_id == "nct01496742_fi" &
        arm %in% c("Placebo + Platinum + Pemetrexed",
                   "MetMAb + Platinum + Pemetrexed") ~ 7L,
      dataset_id == "nct01366131_fi" ~ 8L,
      dataset_id == "nct01493843_fi" &
        arm %in% c("Pictilisib (340 mg) + Carboplatin + Paclitaxel",
                   "Placebo (340 mg) + Carboplatin + Paclitaxel") ~ 9L,
      dataset_id == "nct01493843_fi" &
        arm %in% c("Pictilisib (340 mg) + Carboplatin + Paclitaxel + Bevacizumab",
                   "Placebo (340 mg) + Carboplatin + Paclitaxel + Bevacizumab") ~ 10L,
      dataset_id == "nct01493843_fi" &
        arm %in% c("Pictilisib (260 mg) + Carboplatin + Paclitaxel + Bevacizumab",
                   "Placebo (260 mg) + Carboplatin + Paclitaxel + Bevacizumab") ~ 11L,
      dataset_id == "nct02367781_fi" ~ 12L,
      dataset_id == "nct02367794_fi" ~ 13L,
      dataset_id == "nct02657434_fi" ~ 14L
     )
    )
  
  # Censoring and endpoints
  ## Get maximum followup in RCT by study
  data <- dplyr::left_join(
    data,
    data %>%
      dplyr::filter(source_type == "RCT") %>%
      dplyr::group_by(analysis_num) %>%
      dplyr::summarize(max_os_days = max(os_days)),
    by = "analysis_num"
  )
  
  ## Modify time to event based on maximum followup time
  data <- data %>%
    # Reset os_days to 0 if it is less than 0
    dplyr::mutate(os_days = ifelse(os_days < 0, 0, os_days)) %>%
    
    # Right-censor Flatiron patients at last trial activity
    dplyr::mutate(os_status = ifelse(os_days > max_os_days, 0, os_status)) %>%
    
    # Reset os_days for right-censored flatiron patients 
    dplyr::mutate(os_days = ifelse(os_days >= max_os_days, max_os_days, 
                                   os_days))
  
  # Make some IDs for each separate analysis 
  # (there may be multiple for a given trial/EC combo)
  data <- data %>%
    dplyr::arrange(analysis_num) %>%
    dplyr::group_by(analysis_num, dataset_id) %>%
    dplyr::mutate(dataset_counter = ifelse(dplyr::row_number() == 1L, 1L, 0L)) %>%
    dplyr::group_by(dataset_id) %>%
    dplyr::mutate(
      dataset_counter = cumsum(dataset_counter),
      analysis_id = paste0(dataset_id, "_", dataset_counter),
      dataset_counter = NULL
    ) %>%
    dplyr::ungroup()
  
  # Convert character variables to factors
  data <- data %>% dplyr::mutate(
    sex = factor(sex, levels = c("M", "F")),
    race_grouped = relevel(as.factor(race_grouped), ref = "White"),
    smoker_grouped = relevel(as.factor(smoker_grouped), ref = "Never"),
    stage_grouped = relevel(as.factor(stage_grouped), ref = "Early")
  )
  
  # Return
  return(data)
}

# Summarize analysis sample ----------------------------------------------------
#' New analysis object
#'
#' Create a new analysis object that summarizes all pairwise analyses used for
#' this study. Each pairwise analysis is uniquely identified by the clinical trial 
#' and the experimental vs. comparator arm. 
#' @param data A dataset returned by [ecdata::pin_nsclc()] and preprocessed
#' using [preprocess()].
#' @return An object of class `analysis` that inherits from [`dplyr::tibble`].Each
#' row is a unique analysis. 
#' @export
new_analysis <- function(data) {
  x <- data %>%
    dplyr::group_by(analysis_num, analysis_id, arm_type, arm, source, source_type, 
             dataset_id) %>%
    dplyr::tally() %>%
    tidyr::pivot_wider(names_from = c(source_type, arm_type),
                values_from = c(source, arm, n)) %>%
    dplyr::select(
      -arm_RCT_Comparator,# Same as external control comparator
      -source_RCT_Comparator # Same source as RCT experimental arm
    ) %>% 
    dplyr::rename(source_ec = source_External_Comparator,
                  source_rct = source_RCT_Experimental,
                  arm_comparator = arm_External_Comparator,
                  arm_experimental = arm_RCT_Experimental,
                  n_ec = n_External_Comparator,
                  n_ic = n_RCT_Comparator,
                  n_trt = n_RCT_Experimental) %>%
    dplyr::ungroup()
  class(x) <- c("analysis", class(x))
  return(x)
}

#' HTML analysis table
#' 
#' Create an HTML table using using [knitr::kable()] (and styled with [kableExtra::kable_styling()]) 
#' that summarizes an [`analysis`][new_analysis()] object. For each analysis, it indicates
#' the external control data source, the clinicaltrials.gov identifier for the clinical trial,
#' the experimental and comparator arm, and the sample size for both the external control and
#' the experimental arm of the clinical trial.
#' @export
html_table.analysis <- function(x) {
  x %>%
    dplyr::select(-dataset_id) %>%
    dplyr::arrange(analysis_num) %>%
    knitr::kable(
      col.names = c("Number", "ID", 
                    rep(c("External", "Trial"), 2),
                    "External control", "Internal control", "Experimental arm")
    ) %>%
    kableExtra::kable_styling() %>%
    kableExtra::add_header_above(c(" " = 2, "Source" = 2, "Arm" = 2, "N" = 3))
}

# Inspect missing data ---------------------------------------------------------
#' Count missing observations
#' 
#' Count missing observations for variables contained in a dataset.
#' @inheritParams new_analysis
#' @param vars A character vector of the names of variables to count the 
#' number of missing observations for. If the vector is named, then the names
#' are used as labels for the variables. 
#' 
#' @return An object of class `missing_counts` inheriting from [`dplyr::tibble`] 
#' with the following columns:
#' \describe{
#' \item{label}{The label for each variable in `vars`.}
#' \item{source_type}{Indicates whether the data is from the external control or the trial.}
#' \item{analysis_num}{An integer indicating the analyis number.}
#' \item{n_missing}{The number of missing observations.}
#' \item{n}{The number of total observations.}
#' \item{prop_missing}{The proportion of observations that are missing.}
#' }
#' @export
count_missing <- function(data, vars) {
  lab_lookup <- data.frame(var = vars, label = names(vars))
  vars <- unname(vars)
  
  x <- data %>%
    dplyr::select(one_of(vars, "source_type", "analysis_num")) %>%
    dplyr::mutate_at(vars,
              function (x) ifelse(is.na(x), 1, 0)) %>%
    tidyr::pivot_longer(cols = vars, names_to = "var",
                 values_to = "missing") %>%
    dplyr::left_join(lab_lookup, by = "var") %>%
    dplyr::arrange(analysis_num) %>%
    dplyr::group_by(label, source_type, analysis_num) %>%
    dplyr::summarise(n_missing = sum(missing), n = dplyr::n()) %>%
    dplyr::mutate(prop_missing = n_missing/n) 
  class(x) <- c("missing_counts", class(x))
  return(x)
}

#' Plot missing observations
#' 
#' Bar plot to summarize the number of missing observations for each variable
#' selected in [count_missing()]. Plots are faceted by analysis number and
#' variables are on the x-axis in each facet.    
#' @param object An object of class [`missing_counts`][count_missing()].
#' @param yvar Character vector of length 1 indicating whether the y-axis should 
#' report the proportion of missing observations (`"prop_missing"`) or the
#'  number of missing observations (`"n_missing"`).
#'  
#' @return A [ggplot2::ggplot] object.
#' @export
autoplot.missing_counts <- function(object, yvar = c("prop_missing", "n_missing")) {
  yvar <- match.arg(yvar)
  scales <- if (yvar == "prop_missing") "free_x" else "free"
  y_lab <- if (yvar == "prop_missing") "Proportion missing" else "Number missing"
  
  ggplot(object %>%
           dplyr::filter(prop_missing < 1), 
         aes(x = label, y = .data[[yvar]], fill = source_type)) +
    geom_bar(stat = "identity",  position = "dodge") +
    facet_wrap(~analysis_num, ncol = 2, scales = scales) +
    xlab("") +
    scale_fill_discrete(name = "Source") +
    ylab(y_lab) +
    theme(legend.position = "bottom",
          text = element_text(size = 10))
}

# Distribution of covariates ---------------------------------------------------
#' Density plot for continuous variable
#' 
#' Plot kernel density estimates of a continuous variable faceted by the 
#' analysis number. 
#' @inheritParams count_missing
#' @param var A character vector of  length 1 indicating the variable to plot.
#' @param xlim Limits for the x-axis passed to [ggplot2::coord_cartesian()].
#' 
#' @return A [`ggplot2::ggplot`] object.
#' @export
plot_density.data.frame <- function(data, var, xlim = NULL) {
  pdata <- data %>%
    dplyr::select(dplyr::one_of("source_type", "analysis_num", var)) %>%
    tidyr::pivot_longer(cols = var, names_to = "var", values_to = "value") %>%
    dplyr::filter(!is.na(value))
  
  ggplot(pdata,
         aes(x = value, col = source_type)) +
    geom_density() +
    facet_wrap(~analysis_num, scales = "free") +
    xlab("") + ylab("Density") +
    scale_color_discrete(name = "Source") +
    theme(legend.position = "bottom") +
    coord_cartesian(xlim = xlim)
} 

#' Bar plot for categorical variable
#' 
#' Draw bar plots for a categorical variable faceted by the 
#' analysis number. 
#' @inheritParams plot_density
#' 
#' @return A [`ggplot2::ggplot`] object.
#' @export
plot_bar <- function(data, var) {
  pdata <- data %>%
    dplyr::select(one_of("source_type", "analysis_num", var)) %>%
    tidyr::pivot_longer(cols = var, names_to = "var", values_to = "value") %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::group_by(source_type, analysis_num, var, value) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::group_by(source_type, analysis_num, var) %>%
    dplyr::mutate(prop = n / sum(n))
  
  ggplot(pdata, aes(x = value, y = prop, fill = source_type)) +
    geom_bar(stat = "identity",  position = "dodge") +
    facet_wrap(~analysis_num, scales = "free_x") +
    xlab("") +
    scale_y_continuous(breaks = seq(0, 1, .1)) +
    scale_fill_discrete(name = "Source") +
    ylab("Proportion") +
    theme(legend.position = "bottom")
}

# Small sample sizes -----------------------------------------------------------
#' Counts by group
#' 
#' Count the number of unique values of one ore more variables by groups. 
#' @inheritParams count_missing
#' @param vars A character vector indicating the variables to count unique
#' values for.
#' @param by Columns to group by.
#' @param max_n The maximum number of values in a group in the returned data.
#' 
#' @return A [`dplyr::tibble`].
#' @export
count_by <- function(data, vars, by, max_n = Inf) {
  data %>%
    dplyr::select(dplyr::one_of(by, vars)) %>%
    tidyr::pivot_longer(cols = vars, names_to = "var", values_to = "value") %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::group_by(dplyr::across(c(by, "var", "value"))) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::filter(n <= max_n) %>%
    dplyr::arrange(dplyr::across(by))
}
