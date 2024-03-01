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
#' \item{predicted}{A vector of model prediction values. For continuous outcome, this is a vector of model scores; 
#' for binary outcome, this is a vector representing the probability of the positive class;
#' for time to event outcome, this is a vector of risk scores}
#' @author Aixiang Jiang
#' @import xgboost
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

XGBtraining_predict = function(xgbtrainingObj = NULL, newdata = NULL, newY = FALSE, outfile = "nameWithPath") {
  if(is.null(xgbtrainingObj)){
    stop("XGBtraining is null")
  }
  
  testdat = newdata[,xgbtrainingObj$XGBoost_model$feature_names]
  test = xgboost::xgb.DMatrix(as.matrix(testdat)) ## change to as.matrix from data.matrix to deal with potential one variable situation
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