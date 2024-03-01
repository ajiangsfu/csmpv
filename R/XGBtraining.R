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

XGBtraining = function(data, biomks = NULL,outcomeType = c("binary","continuous","time-to-event"), Y =NULL, time=NULL, event=NULL, nrounds = 5,
                       nthread = 2, gamma = 1, max_depth = 3, eta = 0.3, outfile = "nameWithPath") {
  outcomeType = outcomeType[1]
  modeln = NA
  
  ## add one more on 20230705 for cox
  h0 = NULL ## cumulative baseline hazard table
  
  if(outcomeType == "binary"){
    x.train = data[,c(biomks, Y)]
    ### XGB
    num_feature = dim(x.train)[2]
    x.train.xgb = data.matrix(x.train)
    ## add as.matrix to deal with potential one variable situation
    dtrain = list(data=as.matrix(x.train.xgb[,c(1:(num_feature-1))]),label=x.train.xgb[,num_feature])	
    Dtrain = xgboost::xgb.DMatrix(dtrain$data,label=dtrain$label)
    modeln = xgboost::xgboost(
      objective = "binary:logistic",
      data = Dtrain,
      nrounds = nrounds, # max number of boosting iterations
      nthread = nthread, ## Number of parallel threads used to run XGBoost, this is not important for small data set
      verbose = 2, 
      # If 0, xgboost will stay silent. If 1, xgboost will print information of performance. 
      # If 2, xgboost will print information of both performance and construction progress information
      gamma = gamma, # Minimum loss reduction required to make a further partition on a leaf node of the tree
      max_depth = max_depth, # Maximum depth of a tree, default is 6
      eta = eta # this is default, Step size shrinkage used in update to prevents overfitting
    )
  }else if(outcomeType == "continuous"){
    x.train = data[,c(biomks, Y)]
    ### XGB
    num_feature = dim(x.train)[2]
    x.train.xgb = data.matrix(x.train)
    dtrain = list(data=x.train.xgb[,c(1:(num_feature-1))],label=x.train.xgb[,num_feature])	
    Dtrain = xgboost::xgb.DMatrix(dtrain$data,label=dtrain$label)
    modeln = xgboost::xgboost(
      objective = "reg:squarederror",
      data = Dtrain,
      nrounds = nrounds, # max number of boosting iterations
      nthread = nthread, ## Number of parallel threads used to run XGBoost, this is not important for small data set
      verbose = 2, 
      # If 0, xgboost will stay silent. If 1, xgboost will print information of performance. 
      # If 2, xgboost will print information of both performance and construction progress information
      gamma = gamma, # Minimum loss reduction required to make a further partition on a leaf node of the tree
      max_depth = max_depth, # Maximum depth of a tree, default is 6
      eta = eta # this is default, Step size shrinkage used in update to prevents overfitting
    )
  }else if(outcomeType == "time-to-event"){
    x.train = data[,c(biomks, time, event)]
    ### XGB
    num_feature = dim(x.train)[2]
    x.train.xgb = data.matrix(x.train)
    ## essential, y = -time if event = 0, and y = time if event = 1
    dtrain = list(data=x.train.xgb[,c(1:(num_feature-2))],label=x.train.xgb[,(num_feature-1)]*(-(-1)^(as.numeric(x.train.xgb[,num_feature]))))	
    Dtrain = xgboost::xgb.DMatrix(dtrain$data,label=dtrain$label)
    modeln = xgboost::xgboost(
      objective = "survival:cox",
      data = Dtrain,
      nrounds = nrounds, # max number of boosting iterations
      nthread = nthread, ## Number of parallel threads used to run XGBoost, this is not important for small data set
      verbose = 2, 
      # If 0, xgboost will stay silent. If 1, xgboost will print information of performance. 
      # If 2, xgboost will print information of both performance and construction progress information
      gamma = gamma, # Minimum loss reduction required to make a further partition on a leaf node of the tree
      max_depth = max_depth, # Maximum depth of a tree, default is 6
      eta = eta # this is default, Step size shrinkage used in update to prevents overfitting
    )
    coxfit = survival::coxph(Surv(data[,time], data[, event]) ~ 1)
    h0 = basehaz(coxfit)
  }else{
    stop("Please select the correct outcome type")
  }
  
  ## write out the internal validation results for each boosting iteration
  sink(paste0(outfile,"_Internal_validation.txt"))
  print(modeln$evaluation_log)
  sink()
  
  
  pn = stats::predict(modeln, Dtrain)
  names(pn) = rownames(data)
  
  #outs = list(modeln, pn, data[,Y], outcomeType)
  outs = list(modeln, pn, Y, outcomeType)
  names(outs) = c("XGBoost_model", "XGBoost_score", "Y", "outcomeType")
  
  ## add more on 20230705 and modify on 20230711
  if(!is.null(h0)){
    #outs = list(modeln, pn, h0, data[,time], data[,event], outcomeType)
    outs = list(modeln, pn, h0, time, event, outcomeType)
    names(outs) = c("XGBoost_model", "XGBoost_score", "h0", "time", "event", "outcomeType")
  }
  
  return(outs)
}