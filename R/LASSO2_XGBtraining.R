#' Variable Selection with LASSO2 and Modeling with XGBoost
#' 
#' @description This function performs a two-step process: variable selection using LASSO2 and building a predictive model using XGBoost. 
#' @details
#' The first part of LASSO2_XGBtraining involves variable selection with LASSO2, typically based on the mean lambda.1se from 10 iterations of n-fold 
#' cross-validation-based LASSO regression. In each iteration, a lambda.1se refers to the largest value of lambda such that the error is within
#' 1 standard error of the minimum. However, if only one or no variable is selected, the cross-validation results are ignored, 
#' and a minimum of two remaining variables is ensured through full-data lambda simulations.
#' 
#' The second part of LASSO2_XGBtraining involves ignoring the shrunk LASSO coefficients and building an XGBoost model.
#' It is suitable for three types of outcomes: continuous, binary, and time-to-event.
#' 
#' @param data A data matrix or data frame containing samples in rows and features/traits in columns.
#' @param standardization A logical value indicating if standardization is needed before variable selection. Default is FALSE.
#' @param columnWise A logical value indicating if column-wise or row-wise normalization is needed for standardization. Default is TRUE. 
#' This parameter is only meaningful when standardization is TRUE.
#' @param biomks A vector of potential biomarkers for variable selection. These should be a subset of the column names in the data parameter.
#' @param outcomeType The type of the outcome variable: "binary" (default), "continuous", or "time-to-event".
#' @param Y The name of the outcome variable when the outcome type is either "binary" or "continuous".
#' @param time The name of the time variable when the outcome type is "time-to-event".
#' @param event The name of the event variable when the outcome type is "time-to-event".
#' @param nfolds The number of folds for cross-validation. The default is 10.
#' @param nrounds The maximum number of boosting iterations for the XGBoost model.
#' @param nthread The number of parallel threads used for running XGBoost.
#' @param gamma The minimum loss reduction required to make a further partition on a leaf node of the tree.
#' @param max_depth The maximum depth of a tree in the XGBoost model.
#' @param eta The learning rate for the XGBoost model.
#' @param outfile A string for the output file, including the path if necessary, but without the file type extension.
#' @author Aixiang Jiang
#' @return A list is returned:
#' \item{XGBoost_model}{An XGBoost model}
#' \item{XGBoost_model_score}{Model scores for the given training data set. 
#'  For a continuous outcome variable, this is a vector of the estimated continuous values; 
#'  for a binary outcome variable, this is a vector representing the probability of the positive class;
#'  for time-to-event outcome, this a vector of risk scores}
#' @references 
#'   Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent (2010), Journal of Statistical Software, Vol. 33(1), 1-22, doi:10.18637/jss.v033.i01.
#'   
#'   Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol. 39(5), 1-13, doi:10.18637/jss.v039.i05.
#'   
#'   Tianqi Chen and Carlos Guestrin, "XGBoost: A Scalable Tree Boosting System", 22nd SIGKDD Conference on Knowledge Discovery and Data Mining, 2016, https://arxiv.org/abs/1603.02754
#' @examples
#' # Load in data sets:
#' data("datlist", package = "csmpv")
#' tdat = datlist$training
#' 
#' # The function saves files locally. You can define your own temporary directory. 
#' # If not, tempdir() can be used to get the system's temporary directory.
#' temp_dir = tempdir()
#' # As an example, let's define Xvars, which will be used later:
#' Xvars = c("highIPI", "B.Symptoms", "MYC.IHC", "BCL2.IHC", "CD10.IHC", "BCL6.IHC")

#' # The function can work with three different outcome types. 
#' # Here, we use binary as an example:
#' blxfit = LASSO2_XGBtraining(data = tdat, biomks = Xvars, Y = "DZsig",
#'                            outfile = paste0(temp_dir, "/binary_LASSO2_XGBoost"))
#' # You might save the files to the directory you want.
#' 
#' # To delete the "temp_dir", use the following:
#' unlink(temp_dir)


#' @export

LASSO2_XGBtraining = function(data = NULL, standardization = FALSE, columnWise = TRUE, biomks = NULL, outcomeType = c("binary","continuous","time-to-event"), 
                             Y = NULL, time = NULL, event = NULL, nfolds = 10, nrounds = 5, nthread = 2, gamma = 1, max_depth = 3, eta = 0.3, outfile = "nameWithPath") {
  ## I could call LASSO directly, but to be more confident, I call cv.glmnet instead
  if(is.null(data)){
    stop("Please input a data set")
  }
  if(standardization){
    data = standardize(data, byrow = !columnWise)
  }
  
  outcomeType = outcomeType[1]

  lasres = LASSO2(data = data, biomks = biomks, outcomeType = outcomeType, 
                   Y = Y, time = time, event = event, nfolds = nfolds, outfile = outfile)
  svars = names(lasres$coefs)
  xgbout = NULL
  if(length(svars)>1){## it does not work if only there is no or one variable
    ## call XGBtraining
    xgbout = XGBtraining(data = data, biomks = svars,outcomeType = outcomeType, Y = Y, time = time, event = event, nrounds = nrounds,
                         nthread = nthread, gamma = gamma, max_depth = max_depth, eta = eta, outfile = outfile) 
  }
  ## return output of XGBtraining
  return(xgbout)

}


