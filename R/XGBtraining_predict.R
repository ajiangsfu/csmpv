#' Predicting XGBoost Model Scores and Performing Validation
#' @description 
#' This function predicts XGBoost model scores using an XGBtraining object and a new dataset. 
#' It converts the input data to the required xgb.DMatrix format and returns the model scores. 
#' If the new dataset includes an outcome variable, the function also performs validation, 
#' comparing predictions with observed outcomes.
#' @param xgbtrainingObj An XGBtraining object returned from the XGBtraining function.
#' @param newdata A data matrix or a data frame, samples are in rows, and features/traits are in columns.
#' @param newY A logical variable indicating if 'newdata' contains the outcome variable.
#' @param outfile A string for the output file including path if necessary but without file type extension. 
#' @return A vector of predicted values is return. If an outcome variable is available for the new dataset, validation is performed.
#'         For continuous outcome, this is a vector of model scores. 
#'         For binary outcome, this is a vector representing the probability of the positive class.
#'         For time to event outcome, this is a vector of risk scores.
#' @author Aixiang Jiang
#' @import xgboost
#' @import survival
#' @importFrom stats .getXlevels
#' @importFrom stats model.matrix
#' @importFrom stats terms
#' @importFrom stats model.frame

#' @references 
#'   Tianqi Chen and Carlos Guestrin, "XGBoost: A Scalable Tree Boosting System", 22nd SIGKDD Conference on Knowledge Discovery and Data Mining, 2016, https://arxiv.org/abs/1603.02754
#'   
#'   Harrell Jr F (2023). rms: Regression Modeling Strategies_. R package version 6.7-1, <https://CRAN.R-project.org/package=rms>
#'   
#'   Harrell Jr F (2023). Hmisc: Harrell Miscellaneous_. R package version 5.1-1, <https://CRAN.R-project.org/package=Hmisc>

#' @examples
#' # Load in data sets:
#' data("datlist", package = "csmpv")
#' tdat = datlist$training
#' vdat = datlist$validation
#' 
#' # The function saves files locally. You can define your own temporary directory. 
#' # If not, tempdir() can be used to get the system's temporary directory.
#' temp_dir = tempdir()
#' # As an example, let's define Xvars, which will be used later:
#' Xvars = c("highIPI", "B.Symptoms", "MYC.IHC", "BCL2.IHC", "CD10.IHC", "BCL6.IHC")
#'
#' # The function can work with multiple models and multiple outcome types. 
#' # Here, we provide an example using the XGBoost model with a time-to-event outcome:
#' txfit = XGBtraining(data = tdat, biomks = Xvars,
#'                     outcomeType = "time-to-event",
#'                     time = "FFP..Years.",event = "Code.FFP",
#'                     outfile = paste0(temp_dir, "/survival_XGBoost"))
#' ptxfit = XGBtraining_predict(txfit, newdata = vdat,
#'                     outfile = paste0(temp_dir, "/pred_XGBoost_time_to_event"))
#' # To delete the "temp_dir", use the following:
#' unlink(temp_dir)


#' @export

XGBtraining_predict <- function(xgbtrainingObj = NULL, newdata = NULL, newY = FALSE, outfile = "nameWithPath") {
  
  # ========== INPUT VALIDATION ==========
  if (is.null(xgbtrainingObj)) {
    stop("xgbtrainingObj cannot be NULL")
  }
  if (is.null(newdata)) {
    stop("newdata cannot be NULL")
  }
  if (!inherits(newdata, "data.frame")) {
    stop("newdata must be a data.frame")
  }
  
  # ========== EXTRACT MODEL COMPONENTS ==========
  model <- xgbtrainingObj$XGBoost_model
  outcomeType <- xgbtrainingObj$outcomeType
  
  if (is.null(model)) {
    stop("XGBoost model not found in xgbtrainingObj")
  }
  
  # Extract stored encoding information
  factor_info <- attr(model, "factor_info")
  stored_formula <- attr(model, "formula")
  stored_feature_names <- attr(model, "feature_names")
  
  # Backward compatibility
  if (is.null(stored_feature_names)) {
    stored_feature_names <- model$feature_names
    if (is.null(stored_feature_names)) {
      stored_feature_names <- xgbtrainingObj$feature_names
    }
  }
  
  if (is.null(stored_feature_names)) {
    stop("Feature names not found in model object")
  }
  
  # Create the matrix
  test_matrix <- tryCatch({
    create_consistent_matrix(newdata, stored_formula, factor_info, stored_feature_names)
  }, error = function(e) {
    # Fallback: try with feature names only
    warning("Using fallback method for matrix creation: ", e$message)
    formula_fallback <- as.formula(paste("~", paste(stored_feature_names, collapse = " + ")))
    test_matrix_fallback <- model.matrix(formula_fallback, data = newdata)[, -1, drop = FALSE]
    matrix(as.numeric(test_matrix_fallback), nrow = nrow(test_matrix_fallback))
  })
  
  # ========== PREDICT ==========
  test_dmatrix <- xgboost::xgb.DMatrix(data = test_matrix)
  scores <- stats::predict(model, test_dmatrix)
  names(scores) <- rownames(newdata)
  
  # ========== VALIDATION IF REQUESTED ==========
  if (newY) {
    # Extract outcome information
    Y <- xgbtrainingObj$Y
    time <- xgbtrainingObj$time
    event <- xgbtrainingObj$event
    baseHz <- xgbtrainingObj$h0
    
    # Validate we have the required information
    if (outcomeType %in% c("binary", "continuous") && is.null(Y)) {
      warning("Y not found for validation, skipping validation")
      newY <- FALSE
    }
    if (outcomeType == "time-to-event" && (is.null(time) || is.null(event))) {
      warning("time and/or event not found for validation, skipping validation")
      newY <- FALSE
    }
    
    if (newY) {
      if (outcomeType == "continuous") {
        outs <- validation(predicted = scores, outcomeType = "continuous", 
                           trueY = newdata[[Y]], outfile = outfile)
      } else if (outcomeType == "binary") {
        outs <- validation(predicted = scores, outcomeType = "binary", 
                           trueY = newdata[[Y]], outfile = outfile)
      } else if (outcomeType == "time-to-event") {
        outs <- validation(predicted = scores, outcomeType = "time-to-event", 
                           time = newdata[[time]], trueEvent = newdata[[event]], 
                           baseHz = baseHz, outfile = outfile)
      }
      return(outs)
    }
  }
  
  # Return just scores if no validation requested
  return(scores)
}