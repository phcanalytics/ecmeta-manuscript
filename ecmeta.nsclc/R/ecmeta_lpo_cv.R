#' Meta-analytic leave-p-out cross validation
#' 
#' Use leave-p-out cross-validation (LpO CV) to compare an external control
#' analysis approach that uses the `ecmeta` package to adjust for bias to a 
#' standard approach that makes no adjustments. LpO CV uses p data points in the 
#' data for testing and all other data points for training.
#' 
#' @param data A data frame with columns `param`, `analysis_num`,
#'  `estimate`, and `se`. `param` should have three values representing the log
#'  hazard ratio of a comparison between (i) the treatment and the external 
#'  control (`trt_ec`), (ii) the internal control and the external control
#'  (`ic_ec`), and (iii) the treatment and the internal control (`trt_ic`). 
#'  `analysis_num` is the analysis (i.e., study) number, `estimate` is the point
#'  estimate, and `se` is the standard error. 
#' @param p Number of data points to use for testing. Default is 1 which is
#' leave-one-out cross validation. 
#' @param method The method/software used to estimate the meta-analytic model
#' as in `ecmeta`.
#' @param ecmeta_args A named list containing additional arguments to pass
#' to [ecmeta::ecmeta()].
#' @param predict_args A named list containing additional arguments to pass to
#' [ecmeta::predict.ecmeta()].
#' 
#' @details LpO CV requires training and validating
#' the dataset \eqn{C_p^n} where \eqn{n} is the number of observations in the original dataset.
#' For each of the \eqn{C_p^n} iterations, a meta-analytic model is fit to the training 
#' data using `ecmeta::ecmeta()` and a prediction on the training set is made with 
#' `ecmeta::predict.ecmeta()`. Predictions that are adjusted with the meta-analytic method 
#' are compared to those without using any adjustment.
#' 
#' @return A data frame where there are `np` rows (`n` = each time the model
#' is trained and tested; `p` is the number of testing observations for each
#' of the `n` iterations). Each row contains performance metrics.
#' @export
ecmeta_lpo_cv <- function(data, p = 1, method = "ml", ecmeta_args = NULL,
                          predict_args = NULL) {
  
  if (!"ecmeta" %in% rownames(installed.packages())) {
    stop("You must install 'ecmeta'.")
  }
  
  id <- unique(data$analysis_num)
  test_ids <- utils::combn(id, m = p, simplify = FALSE)
  metrics <- lapply(test_ids, ecmeta_lpo_cv1, data = data, method = method, 
                    ecmeta_args = ecmeta_args, predict_args = predict_args)
  metrics <- dplyr::bind_rows(metrics, .id = "iteration")
  metrics
}


ecmeta_lpo_cv1 <- function(test_ids, data, method, ecmeta_args,
                           predict_args) {
  
  # Fit meta-analytic model on training set
  train <- data[!data$analysis_num %in% test_ids, ]
  loghr_ic_ec <- ecmeta::as_loghr_data(train[train$param == "ic_ec", ],
                                       estimate = "estimate",
                                       standard_error = "se")
  loghr_ecmeta <- do.call(
    ecmeta::ecmeta,
    args = c(list(data = loghr_ic_ec,
                  method = method),
             ecmeta_args)
  )
  
  # Predictions on the test set
  n_test_ids <- length(test_ids)
  y <- yhat_adj <- yhat_noadj <- bias <-
    residuals_adj <- std_residuals_adj <- 
    residuals_noadj <- std_residuals_noadj <- rep(NA, n_test_ids)
  
  for (j in 1:n_test_ids) {
    ## Predictions on test set
    test <- data[data$analysis_num == test_ids[j], ]
    new_loghr_trt_ec <- ecmeta::as_loghr_data(test[test$param == "trt_ec", ],
                                              estimate = "estimate",
                                              standard_error = "se")
    loghr_new <- do.call(
      predict,
      args = c(list(object = loghr_ecmeta,
                    newdata = new_loghr_trt_ec),
               predict_args)
    )
    
    ## Compute metrics
    new_loghr_trt_ic <- test[test$param == "trt_ic", ]
    y[j] <- new_loghr_trt_ic$estimate
    bias[j] <- mean(loghr_new$loghr[, "trt_ic"] - loghr_new$loghr[, "trt_ec"])
    
    ### Using meta-analytic method to adjust
    residuals_adj_post <- new_loghr_trt_ic$estimate - loghr_new$loghr[, "trt_ic"]
    residuals_adj_sd <- sqrt(var(loghr_new$loghr[, "trt_ic"]))
    yhat_adj[j] <-  median(loghr_new$loghr[, "trt_ic"])
    residuals_adj[j] <- mean(residuals_adj_post)
    std_residuals_adj[j] <- residuals_adj[j]/residuals_adj_sd
    
    ### Standard approach with no adjustment
    residuals_noadj_post <- new_loghr_trt_ic$estimate - loghr_new$loghr[, "trt_ec"]
    residuals_noadj_sd <- sqrt(var(loghr_new$loghr[, "trt_ec"]))
    yhat_noadj[j] <-  median(loghr_new$loghr[, "trt_ec"])
    residuals_noadj[j] <- mean(residuals_noadj_post)
    std_residuals_noadj[j] <- residuals_noadj[j]/residuals_noadj_sd
  }
  
  # Return outcomes
  out <- data.frame(test_id = rep(test_ids, times = 2),
                    adjust = rep(c("Yes", "No"), each = n_test_ids),
                    y = rep(y, times = 2),
                    bias = rep(bias, times = 2),
                    yhat = c(yhat_adj, yhat_noadj),
                    residual = c(residuals_adj, residuals_noadj),
                    std_residual = c(std_residuals_adj, std_residuals_noadj))
}

