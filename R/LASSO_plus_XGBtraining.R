#' LASSO_plus_XGBtraining: Variable Selection and XGBoost Modeling
#' 
#' @description
#' This function performs variable selection using LASSO_plus, followed by building an XGBoost model.
#' LASSO_plus is a method for variable selection, and XGBoost is a gradient boosting algorithm for modeling.
#'
#' @param data A data matrix or a data frame where samples are in rows, and features/traits are in columns.
#' @param standardization A logical variable to indicate if standardization is needed before variable selection. The default is FALSE.
#' @param columnWise A logical variable to indicate if column-wise or row-wise normalization is needed. The default is TRUE, 
#' which is to do column-wise normalization. This is only meaningful when "standardization" is TRUE.
#' @param biomks A vector of potential biomarkers for variable selection. They should be a subset of "data" column names. 
#' @param outcomeType Outcome variable type. There are three choices: "binary" (default), "continuous", and "time-to-event".  
#' @param Y Outcome variable name when the outcome type is either "binary" or "continuous". 
#' @param time Time variable name when outcome type is "time-to-event".
#' @param event Event variable name when outcome type is "time-to-event".
#' @param topN An integer to indicate how many variables we intend to select.
#' @param nrounds Max number of boosting iterations.
#' @param nthread Number of parallel threads used to run XGBoost.
#' @param gamma Minimum loss reduction required to make a further partition on a leaf node of the tree.
#' @param max_depth Maximum depth of a tree. Increasing this value will make the model more complex and more likely to overfit.
#' @param outfile A string representing the output file, including the path if necessary, but without the file type extension.
#' @param eta The learning rate for the XGBoost model.
#' @param height An integer to indicate the forest plot height in inches.
#' 
#' @return A list is returned: 
#' \item{XGBoost_model}{An XGBoost model}
#' \item{XGBoost_model_score}{Model scores for the given training data set. 
#'  For a continuous outcome variable, this is a vector of the estimated continuous values; 
#'  for a binary outcome variable, this is a vector representing the probability of the positive class;
#'  for time-to-event outcome, this is a vector of risk scores}
#' @author Aixiang Jiang

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
#' # Here, we use continuous as an example:
#' clpxfit = LASSO_plus_XGBtraining(data = tdat, biomks = Xvars,
#'           outcomeType = "continuous", Y = "Age", topN = 5,
#'           outfile = paste0(temp_dir, "/continuous_LASSO_plus_XGBoost"))
#' # You might save the files to the directory you want.
#' 
#' # To delete the "temp_dir", use the following:
#' unlink(temp_dir)

#' @export

LASSO_plus_XGBtraining = function(data = NULL, standardization = FALSE, columnWise = TRUE, biomks = NULL, outcomeType = c("binary","continuous","time-to-event"), 
                             Y = NULL, time = NULL, event = NULL, topN = 10, nrounds = 5, nthread = 2, gamma = 1, max_depth = 3, eta = 0.3, 
                             outfile = "nameWithPath",height = 6) {
  
  ## call LASSO_plus
  afit = LASSO_plus(data = data, standardization = standardization, columnWise = columnWise, biomks = biomks, 
                 outcomeType = outcomeType, Y = Y, time = time, event = event, topN = topN, outfile = outfile, height = height)
  
  ## get updated biomks
  aform = afit$fit$formula
  tmp = strsplit(toString(aform), split = ",")[[1]]
  tmp = tmp[length(tmp)]
  tmp = strsplit(tmp, split = "\\+")[[1]]
  if(length(tmp) > 1){
    tmp = sapply(tmp, function(xx){
      gsub(" ", "", xx)
    })
  }else{
    tmp = gsub(" ", "", tmp)
  }
  biomks = tmp
  
  ## call XGBtraining
  xgbout = XGBtraining(data = data, biomks = biomks,outcomeType = outcomeType, Y =Y, time = time, event = event, nrounds = nrounds,
                         nthread = nthread, gamma = gamma, max_depth = max_depth, eta = eta, outfile = outfile) 
  
  ## return output of XGBtraining
  return(xgbout)

}


