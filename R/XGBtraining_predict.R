#' Predicting XGBoost Model Scores and Performing Validation
#' 
#' @description This function is designed for predicting XGBoost model scores based on an xgbtrainingObj object and a new dataset. It converts the provided new data format into the required xgb.DMatrix format and returns the corresponding model scores.
#' 
#' If the new dataset contains an outcome variable, the function also performs a validation step, comparing the predictions with the observed outcomes.
#'
#' @param xgbtrainingObj A xgbtrainingObj object
#' @param newdata A data matrix or a data frame, samples are in rows, and features/traits are in columns.
#' @param newY A logical variable indicating if 'newdata' contains the outcome variable.
#' @param outfile A string for the output file including path if necessary but without file type extension. 
#' @return A vector of predicted values is return. If outcome variable is variable for the new data set, validation is processed and
#'  a list with the following items is returned:
#' \item{predicted}{A vector of model prediction values. For continuous outcome, this is a vector of model scores; 
#' for binary outcome, this is a vector representing the probability of the positive class;
#' for time to event outcome, this is a vector of risk scores}
#' @author Aixiang Jiang
#' @import xgboost
#' @references 
#'   Tianqi Chen and Carlos Guestrin, "XGBoost: A Scalable Tree Boosting System", 22nd SIGKDD Conference on Knowledge Discovery and Data Mining, 2016, https://arxiv.org/abs/1603.02754
#'   Harrell Jr F (2023). rms: Regression Modeling Strategies_. R package version 6.7-1, <https://CRAN.R-project.org/package=rms>
#'   Harrell Jr F (2023). Hmisc: Harrell Miscellaneous_. R package version 5.1-1, <https://CRAN.R-project.org/package=Hmisc>

#' @export

XGBtraining_predict = function(xgbtrainingObj = NULL, newdata = NULL, newY = FALSE, outfile = "nameWithPath") {
  if(is.null(xgbtrainingObj)){
    stop("XGBtraining is null")
  }
  testdat = newdata[,xgbtrainingObj$XGBoost_model$feature_names]
  test = xgboost::xgb.DMatrix(data.matrix(testdat))
  scores = stats::predict(xgbtrainingObj$XGBoost_model, test) 
  ## use default, for continuous, this is model score; for binary, this is prob of the positive class; 
  ##   for time to event, this is risk score
  names(scores) = rownames(newdata)
  
  baseHz = xgbtrainingObj$h0
  Y = xgbtrainingObj$Y
  time = xgbtrainingObj$time
  event = xgbtrainingObj$event
  outcomeType = xgbtrainingObj$outcomeType
  
  if(is.null(newdata)){
    stop("Please input a data set")
  }
 
  outs = list(scores)
  
  ## make some basic transformation and call validation function for each outcome type
  if(outcomeType == "continuous"){
    ## validation step: only when newY = true, otherwise, return: outs = list(model_score)
    if(newY){
      outs = validation(predicted = scores, outcomeType = "continuous", trueY = newdata[,Y], outfile = outfile)
    }
    
  }else if(outcomeType == "binary"){
    if(newY){
      outs = validation(predicted = scores, outcomeType = "binary", trueY = newdata[,Y], outfile = outfile)
    }
  }else if(outcomeType == "time-to-event"){
    if(newY){
      outs = validation(predicted = scores, outcomeType = "time-to-event", time = newdata[,time], trueEvent = newdata[,event], 
                        baseHz = baseHz, outfile = outfile)
    }
  }
  return(outs) 

}