# Simulate survival data -------------------------------------------------------

# Function to generate random data:
# We use the non-parametric property of cox models with respect to time to allow 
# us to use exponential distributions
# The rates for the 3 arms are specified by their median time, with coefficients 
# of variation for how they vary between simulated studies
# No censoring is applied. It is assumed all subjects have an event recorded. 
# As the number of events drives the variability of estimates it is not anticipated 
# this would by influential, but should be tested later.
# The number of events can also vary between studies through a coefficient of 
# variation, with min_ne giving a lower limit of the number of events
# A randomisation ratio can be specified which specifies the number of events in 
# the treatment arm to be a multiple of the number of events in the internal 
# control arm in every study. This overrides the ne_trt and cvne_trt options.
# A true hazard ratio can be specified for the treatment arm compared to the 
# internal control in every study. This overrides the med_trt and cv_trt options  
sim_survdata <- function(nstudies = 100, 
                         med_trt = 24, med_ec = 12, med_ic = 15,
                         cv_trt = 0.2, cv_ic = 0.1, cv_ec = 0.3,
                         true_hr = NA,
                         ne_trt = 100, ne_ec = 50, ne_ic = 70, 
                         min_ne = 40,
                         randomisation_ratio = NA, 
                         cvne_trt = 0, cvne_ec = 0, cvne_ic = 0,
                         seed = NULL) {
  
  if(!is.null(seed)) set.seed(seed)
  
  # Dataset where the unit of observation is the study and arm 
  # i.e., there are three rows for each study (IC, EC, and TRT)
  studyarms <- expand.grid(study = 1:nstudies, 
                           med_trt = med_trt, med_ec = med_ec, med_ic=med_ic,
                           cv_trt = cv_trt, cv_ic = cv_ic, cv_ec = cv_ec,
                           ne_trt = ne_trt, ne_ec = ne_ec, ne_ic = ne_ic, 
                           cvne_trt = cvne_trt, cvne_ec = cvne_ec, cvne_ic = cvne_ic) %>%
    pivot_longer(-study) %>%
    mutate(stat = word(name, 1, sep = "_")) %>%
    mutate(arm = word(name, 2, sep = "_")) %>%
    select(-name) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(
      med_study = exp(log(med) + cv * rnorm(nrow(.))),
      ne_study = pmax(min_ne, round(exp(log(ne) + cvne * rnorm(nrow(.))))) 
    ) 
  
  if(!is.na(randomisation_ratio)) {
    studyarms <- merge(
        studyarms,
        filter(studyarms, arm == "ic") %>% 
          select(study, ne_study) %>% 
          rename(ne_ic_study = ne_study), 
        by = "study"
      ) %>%
      mutate(ne_study = if_else(arm =="trt", randomisation_ratio * ne_ic_study, ne_study)) %>%
      select(-ne_ic_study)
  }
  
  if(!is.na(true_hr)) {
    studyarms <- merge(
      studyarms,
      filter(studyarms, arm =="ic") %>% 
        select(study, med_study) %>% 
        rename(med_ic_study = med_study), by = "study") %>%
      mutate(med_study = if_else(arm == "trt", med_ic_study/true_hr, med_study)) %>%
      select(-med_ic_study)
  }
  
  # Dataset where study is the unit of observation
  studies <- studyarms %>%
    pivot_longer(-c(study, arm)) %>%
    mutate(longname = paste(name, arm, sep="_")) %>%
    pivot_wider(id_cols = study, names_from = longname, values_from = value) %>%
    mutate(
      true_loghr_trt_ic = log(med_study_ic) - log(med_study_trt),
      true_loghr_trt_ec = log(med_study_ec) - log(med_study_trt),
      true_loghr_ic_ec = log(med_study_ec) - log(med_study_ic)
    )
  
  # Dataset of simulated survival times where the patient is the unit of observation
  # The number of patients in each study is equal to the number of events
  dat <- studyarms[rep(seq_len(nrow(studyarms)), times = studyarms$ne_study), ] %>%
    mutate(
      y = rexp(nrow(.), 1/med_study),
      event = 1,
      factor_arm = factor(arm , levels = c("ec", "ic", "trt"))
    ) 
  
  list(studies = studies, studyarms = studyarms, dat = dat)
}

## Calculate the three pairwise hazard ratios for a study
estimate_hrs1 <- function(.dat) {
  loghr_trt_ic <- coxph(Surv(y, event) ~ relevel(factor_arm, "ic"), 
                        subset = (arm %in% c("trt", "ic")), data = .dat) %>% 
    coef() %>% 
    na.omit()
  loghr_trt_ec <- coxph(Surv(y,event) ~ relevel(factor_arm, "ec"), 
                        subset = (arm %in% c("trt", "ec")), data = .dat) %>% 
    coef() %>% 
    na.omit()
  loghr_ic_ec <- coxph(Surv(y,event) ~ relevel(factor_arm, "ec"), 
                       subset = (arm %in% c("ic", "ec")), data = .dat) %>% 
    coef() %>% 
    na.omit()
  
  data.frame(loghr_trt_ic = loghr_trt_ic, loghr_trt_ec = loghr_trt_ec, 
             loghr_ic_ec=loghr_ic_ec)
}

# Get the hazard ratios for a set of studies
estimate_hrs <- function(.sim) {
  .sim$dat %>%
    split( . , .$study) %>%
    lapply(estimate_hrs1) %>%
    bind_rows() %>%
    mutate(study = 1:nrow(.)) %>%
    merge(.sim$studies)
}

# Function to plot the simulation scenarios
plot_scenarios <- function(.res , titl="") {
  
  nes <- .res %>%
    select(starts_with("ne_study")) %>%
    pivot_longer(cols=starts_with("ne_study")) %>%
    mutate(shortname=word(name,3,3,"_")) %>%
    ggplot() + 
    geom_histogram(aes(x=value , colour=name , fill=name) , bins=100) +
    facet_wrap(~shortname , nrow=2) + 
    xlab("Number Of Events") +
    theme(legend.position = "none") + 
    ggtitle(titl)
  
  hrs <- .res %>% 
    select(study,contains("loghr")) %>%
    pivot_longer(cols=contains("loghr")) %>%
    mutate(arm=word(name,-2,-1,"_")) %>%
    mutate(type=word(name,1,1,"_")) %>%
    pivot_wider(id_cols=c("study","arm") , names_from="type") %>%
    ggplot() + 
    geom_point(aes(x=true , y=loghr , colour=arm)) + 
    facet_wrap(~arm) + 
    xlab("True Log HR") + 
    ylab("Observed Log HR") +
    theme(legend.position = "none") + 
    ggtitle(titl)
  
  list(nes=nes , hrs=hrs)
}

# Run the simulation -----------------------------------------------------------
# From a set of studies apply the methodology 
run_sim <- function(.res , n_ref = 9) {
  
  gc() %>% invisible()
  
  # Split into simulations - sets of reference studies and a test study
  .res <- split_train_test(.res , n_ref = n_ref)
  
  # Then apply methodology
  out_ml <- apply_ecmeta(.res, method = "ml")
  out_invgamma <- apply_ecmeta(.res, method = "jags", prior = "Inverse gamma")
  out_cauchy <- apply_ecmeta(.res, method = "jags", prior = "Half-Cauchy")
  out_unif <- apply_ecmeta(.res, method = "jags", prior = "Uniform")
  out <- rbind(out_ml, out_cauchy, out_invgamma)
  out$n_ref <- n_ref
  out
}

# Split a large set of simulated studies into sets of n_ref training studies 
# and a test study. Excess studies are dropped
split_train_test <- function(.res , n_ref = 9) {
  n_ref <- max(n_ref)
  n_data <- n_ref + 1 
  n_sims <- floor(nrow(.res)/n_data)
  .res <- .res %>%
    slice_head(n = n_sims * n_data) %>%
    mutate(
      type = rep(c(rep("TRAIN", n_ref), "TEST"), n_sims),
      sim = rep(1:n_sims, each = n_data)
    )
}

# Apply adjustment via ecmeta() for one study
apply_ecmeta1 <- function(data, method, prior_scale) {
  
  # (1) Estimate meta-analytic model
  ## Data from reference studies
  ref_studies <- data[data$type == "TRAIN", c("loghr_ic_ec", "ne_study_ic",
                                             "ne_study_ec")]
  ref_data <- loghr_data(
    estimate = ref_studies$loghr_ic_ec,
    standard_error = sqrt(1/ref_studies$ne_study_ic + 1/ref_studies$ne_study_ec)
  )
  
  ## Model estimation
  if (method == "ml") {
    fit <- ecmeta(
      data = ref_data,
      method = "ml",
      hessian = FALSE
    )
  } else if (method == "jags") {
    prior_scale_spec <- switch( # Specify parameters of prior
      prior_scale,
      "Half-Cauchy" = student_t(0, 10, 1),
      "Inverse gamma" = invgamma(0.001, 0.001),
      "Uniform" = uniform(0, 100)
    )
    fit <- ecmeta(
      data = ref_data,
      method = "jags",
      prior_scale = prior_scale_spec,
      quiet = TRUE
    )
  }
  
  # (2) Prediction for a new study
  ## New data
  new_study <- data %>%
    filter(type == "TEST") %>%
    select(loghr_trt_ec, ne_study_ec, ne_study_trt, true_loghr_trt_ic)
  new_data <- loghr_data(
    estimate = new_study$loghr_trt_ec,
    standard_error = sqrt(1/new_study$ne_study_ec + 1/new_study$ne_study_trt)
  )
  pred <- predict(fit, newdata = new_data, quiet = TRUE)
  
  # (3) Operating characteristics
  sp <- summary(pred)
  est <- sp[sp$param == "loghr_trt_ic", "50%"]
  sd <- sp[sp$param == "loghr_trt_ic", "sd"]
  lower95 <- sp[sp$param == "loghr_trt_ic", "2.5%"]
  upper95 <- sp[sp$param == "loghr_trt_ic", "97.5%"]
  truth <- new_study$true_loghr_trt_ic
  bias <- est - truth
  covered95 <- 1 * (truth >= lower95 & truth <= upper95)
  sig <- 1 * (upper95 < 0)
  data.frame(
    est, truth, sd, bias, covered95, sig
  )
}

apply_ecmeta <- function(.res, method, prior_scale = c("None", "Half-Cauchy", "Inverse gamma", "Uniform")) {
  
  prior_scale <- match.arg(prior_scale)
  
  .res %>% 
    split( . , .$sim) %>%
    lapply(apply_ecmeta1, method = method, prior_scale = prior_scale) %>%
    bind_rows() %>%
    mutate(
      bayesian = ifelse(method == "ml", FALSE, TRUE),
      prior = prior_scale,
      method = case_when(
        prior == "None" ~ "Maximum likelihood",
        prior == "Half-Cauchy" ~ "Bayesian: half-Cauchy",
        prior == "Inverse gamma" ~ "Bayesian: inverse gamma",
        prior == "Uniform" ~ "Bayesian: uniform"
      ),
      sim = 1:nrow(.)
    )
}