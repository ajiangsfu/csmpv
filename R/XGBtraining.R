#' A Wrapper Function for xgboost::xgboost

#' @description 
#' This wrapper function streamlines the process of utilizing the xgboost package for model training. 
#' It takes care of converting the data format to xgb.DMatrix, handling xgboost's specific settings, 
#' and invoking xgboost::xgboost. The function is suitable for all three outcome types: binary, 
#' continuous, and time-to-event. It returns both the trained model and the model scores for the 
#' training dataset.
#' 
#' It's important to note that all independent variables (X variables) should already be selected 
#' and in numeric format when passed to this function. Additionally, this function does not perform 
#' variable selection or automatically convert categorical variables to numeric format.

#' @param data A data matrix or a data frame where samples are in rows and features/traits are in columns.
#' @param biomks A vector of potential biomarkers for variable selection. They should be a subset of the column names in the "data" variable.
#' @param outcomeType The outcome variable type. There are three choices: "binary" (default), "continuous", and "time-to-event".
#' @param Y The outcome variable name when the outcome type is either "binary" or "continuous". When Y is binary, it should be in 0-1 format.
#' @param time The time variable name when the outcome type is "time-to-event".
#' @param event The event variable name when the outcome type is "time-to-event".
#' @param nrounds The maximum number of boosting iterations.
#' @param nthread The number of parallel threads used to run XGBoost.
#' @param gamma The minimum loss reduction required to make a further partition on a leaf node of the tree.
#' @param max_depth The maximum depth of a tree. Increasing this value will make the model more complex and more likely to overfit.
#' @param outfile A string for the output file, including the path if necessary but without the file type extension.
#' @param eta The step size shrinkage used in the update to prevent overfitting.

#' @return A list is returned:
#' \item{XGBoost_model}{An XGBoost model}
#' \item{XGBoost_score}{Scores for the given training data set.
#' For a continuous outcome variable, this is a vector of the estimated continuous values;
#' for a binary outcome variable, this is a vector representing the probability of the positive class;
#' for a time-to-event outcome, this is a vector of risk scores.}
#' \item{h0}{Cumulative baseline hazard table, for time to event outcome only.}
#' \item{Y}{The outcome variable name when the outcome type is either "binary" or "continuous".}
#' \item{time}{The time variable name when the outcome type is "time-to-event".}
#' \item{event}{The event variable name when the outcome type is "time-to-event".}

#' @author Aixiang Jiang
#' @import xgboost
#' @import survival
#' @importFrom stats .getXlevels
#' @importFrom stats model.matrix
#' @importFrom stats terms
#' @importFrom stats model.frame
#' 
#' @references 
#' Tianqi Chen and Carlos Guestrin, "XGBoost: A Scalable Tree Boosting System", 22nd SIGKDD Conference on Knowledge Discovery and Data Mining, 2016, https://arxiv.org/abs/1603.02754

#'@examples
#' # Load in data sets:
#' data("datlist", package = "csmpv")
#' tdat = datlist$training

#' # The function saves files locally. You can define your own temporary directory. 
#' # If not, tempdir() can be used to get the system's temporary directory.
#' temp_dir = tempdir()
#' # As an example, let's define Xvars, which will be used later:
#' Xvars = c("highIPI", "B.Symptoms", "MYC.IHC", "BCL2.IHC", "CD10.IHC", "BCL6.IHC")
#'
#' # The function can work with three outcome types. 
#' # Here, we use time-to-event outcome as an example:
#' txfit = XGBtraining(data = tdat, biomks = Xvars,
#'                     outcomeType = "time-to-event",
#'                     time = "FFP..Years.",event = "Code.FFP",
#'                     outfile = paste0(temp_dir, "/survival_XGBoost"))
#' # To delete the "temp_dir", use the following:
#' unlink(temp_dir)


#' @export

XGBtraining = function(data, biomks = NULL, outcomeType = c("binary", "continuous", "time-to-event"), 
                       Y = NULL, time = NULL, event = NULL, nrounds = 5,
                       nthread = 2, gamma = 1, max_depth = 3, eta = 0.3, outfile = "nameWithPath") {
  
  outcomeType = outcomeType[1]
  modeln = NA
  h0 = NULL ## cumulative baseline hazard table
  
  # Input validation
  if (is.null(biomks) || length(biomks) == 0) {
    stop("biomks cannot be NULL or empty")
  }

  if (outcomeType == "binary") {
    # Validation
    if (is.null(Y)) stop("Y (outcome) must be specified for binary outcome")
    if (!Y %in% colnames(data)) stop(paste("Y column", Y, "not found in data"))
    
    # Prepare data
    x.train <- data[, c(biomks, Y), drop = FALSE]
    
    # Create model matrix
    mm_result <- create_model_matrix(x.train, biomks)
    tdat <- mm_result$matrix
    
    # Extract outcome
    outcome_vals <- x.train[[Y]]
    if (!all(outcome_vals %in% c(0, 1))) {
      warning("Binary outcome should be 0/1. Converting automatically.")
      outcome_vals <- as.numeric(as.factor(outcome_vals)) - 1
    }
    
    # Create DMatrix
    Dtrain <- xgboost::xgb.DMatrix(data = tdat, label = outcome_vals)
    
    # Model parameters
    params <- list(
      objective = "binary:logistic",
      nthread = nthread,
      gamma = gamma,
      max_depth = max_depth,
      eta = eta
    )
    
    # Train model
    modeln <- xgboost::xgb.train(
      params = params,
      data = Dtrain,
      nrounds = nrounds,
      verbose = 2
    )
    
    modeln <- structure(modeln, class = c("xgb.Booster", class(modeln)))
    attr(modeln, "factor_info") <- mm_result$factor_info
    attr(modeln, "feature_names") <- colnames(tdat)
    attr(modeln, "formula") <- mm_result$formula

  } else if (outcomeType == "continuous") {
    # Validation
    if (is.null(Y)) stop("Y (outcome) must be specified for continuous outcome")
    if (!Y %in% colnames(data)) stop(paste("Y column", Y, "not found in data"))
    
    # Prepare data
    x.train <- data[, c(biomks, Y), drop = FALSE]
    
    # Create model matrix
    mm_result <- create_model_matrix(x.train, biomks)
    tdat <- mm_result$matrix
    
    # Extract outcome
    outcome_vals <- as.numeric(x.train[[Y]])
    
    # Create DMatrix
    Dtrain <- xgboost::xgb.DMatrix(data = tdat, label = outcome_vals)
    
    # Model parameters
    params <- list(
      objective = "reg:squarederror",
      nthread = nthread,
      gamma = gamma,
      max_depth = max_depth,
      eta = eta
    )
    
    # Train model
    modeln <- xgboost::xgb.train(
      params = params,
      data = Dtrain,
      nrounds = nrounds,
      verbose = 2
    )
    
    modeln <- structure(modeln, class = c("xgb.Booster", class(modeln)))
    attr(modeln, "factor_info") <- mm_result$factor_info
    attr(modeln, "feature_names") <- colnames(tdat)
    attr(modeln, "formula") <- mm_result$formula
    
  } else if (outcomeType == "time-to-event") {
    # Validation
    if (is.null(time) || is.null(event)) {
      stop("Both time and event must be specified for time-to-event outcome")
    }
    if (!all(c(time, event) %in% colnames(data))) {
      missing <- setdiff(c(time, event), colnames(data))
      stop(paste("Missing columns:", paste(missing, collapse = ", ")))
    }
    
    # Prepare data
    x.train <- data[, c(biomks, time, event), drop = FALSE]
    
    # Create model matrix
    mm_result <- create_model_matrix(x.train, biomks)
    tdat <- mm_result$matrix
    
    # Extract time and event
    time_vals <- as.numeric(x.train[[time]])
    event_vals <- as.numeric(x.train[[event]])
    
    # Validate event values
    if (!all(event_vals %in% c(0, 1))) {
      stop("Event column must contain only 0 (censored) and 1 (event) values")
    }
    
    # Create survival label: y = time for events, y = -time for censored
    survival_label <- time_vals * (2 * event_vals - 1)
    
    # Create DMatrix
    Dtrain <- xgboost::xgb.DMatrix(data = tdat, label = survival_label)
    
    # Model parameters
    params <- list(
      objective = "survival:cox",
      nthread = nthread,
      gamma = gamma,
      max_depth = max_depth,
      eta = eta
    )
    
    # Train model
    modeln <- xgboost::xgb.train(
      params = params,
      data = Dtrain,
      nrounds = nrounds,
      verbose = 2
    )
    
    modeln <- structure(modeln, class = c("xgb.Booster", class(modeln)))
    attr(modeln, "factor_info") <- mm_result$factor_info
    attr(modeln, "feature_names") <- colnames(tdat)
    attr(modeln, "formula") <- mm_result$formula
    
    # Calculate baseline hazard
    coxfit <- survival::coxph(survival::Surv(data[[time]], data[[event]]) ~ 1)
    h0 <- survival::basehaz(coxfit)
    
  } else {
    stop("Please select the correct outcome type: binary, continuous, or time-to-event")
  }
  
  ## Write out the internal validation results for each boosting iteration
  if (!is.null(outfile) && outfile != "") {
    sink(paste0(outfile, "_Internal_validation.txt"))
    print(modeln$evaluation_log)
    sink()
  }
  
  ## Predict on training data
  pn <- stats::predict(modeln, Dtrain)
  names(pn) <- rownames(data)
  
  ## Prepare output list
  outs <- list(
    XGBoost_model = modeln,
    XGBoost_score = pn,
    Y = Y,
    outcomeType = outcomeType
  )
  
  ## Add survival-specific outputs
  if (outcomeType == "time-to-event") {
    outs$h0 <- h0
    outs$time <- time
    outs$event <- event
  }
  
  return(outs)
}