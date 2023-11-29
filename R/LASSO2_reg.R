#' LASSO2 Variable Selection and Regular Regression Modeling
#' 
#' @description
#' This function performs variable selection with LASSO2 but ignores the shrunk LASSO coefficients and builds a regular regression model.
#' 
#' @details
#' The first part of LASSO2_reg involves variable selection with LASSO2, which is typically performed using cross-validation-based LASSO regression. 
#' However, when only one or no variable is selected, the cross-validation results are ignored, and a minimum of two remaining variables 
#' is ensured through full-data lambda simulations.
#' The second part of LASSO2_reg involves ignoring the shrunk LASSO coefficients and building a regular regression model.
#' It is suitable for three types of outcomes: continuous, binary, and time-to-event.
#' @param data A data matrix or a data frame where samples are in rows, and features/traits are in columns.
#' @param standardization A logical variable to indicate if standardization is needed before variable selection. The default is FALSE.
#' @param columnWise A logical variable to indicate if column-wise or row-wise normalization is needed. The default is TRUE, which performs column-wise normalization. 
#'        This is only meaningful when "standardization" is TRUE.
#' @param biomks A vector of potential biomarkers for variable selection. They should be a subset of column names in "data". 
#' @param outcomeType The outcome variable type. There are three choices: "binary" (default), "continuous", and "time-to-event".
#' @param Y The outcome variable name when the outcome type is either "binary" or "continuous".
#' @param time The time variable name when the outcome type is "time-to-event".
#' @param event The event variable name when the outcome type is "time-to-event".
#' @param nfolds The number of folds for cross-validation. The default value is 10.
#' @param outfile A string for the output file including the path, if necessary, but without the file type extension.
#' @importFrom glmnet cv.glmnet
#' @importFrom Matrix nnzero
#' @author Aixiang Jiang
#' 
#' @return 
#' A list is returned with the same output as from confirmVars.
#' \item{fit}{A model with selected variables for the given outcome variable.}
#' \item{allplot}{A list with all plots.}
#' 
#' @references 
#' Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent (2010), Journal of Statistical Software, Vol. 33(1), 1-22, doi:10.18637/jss.v033.i01.
#' Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol. 39(5), 1-13, doi:10.18637/jss.v039.i05.
#' Hastie, T. J. and Pregibon, D. (1992) Generalized linear models. Chapter 6 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
#' Therneau, T., Grambsch, P., Modeling Survival Data: Extending the Cox Model. Springer-Verlag, 2000.
#' Kassambara A, Kosinski M, Biecek P (2021). survminer: Drawing Survival Curves using 'ggplot2'_. R package version 0.4.9,
#' @export

LASSO2_reg = function(data = NULL, standardization = FALSE, columnWise = TRUE, biomks = NULL, outcomeType = c("binary","continuous","time-to-event"), 
                      Y = NULL, time = NULL, event = NULL, nfolds = 10, outfile = "nameWithPath"){
  
  if(is.null(data)){
    stop("Please input a data set")
  }
  if(standardization){
    data = standardize(data, byrow = !columnWise)
  }
  
  alls = NA
  
  outcomeType = outcomeType[1]
  
  lasres = LASSO2(data = data, biomks = biomks, outcomeType = outcomeType, 
                 Y = Y, time = time, event = event, nfolds = nfolds, outfile = outfile)
  svars = names(lasres$coefs)
  
  if(length(svars) > 0){
    fit = confirmVars(data = data, standardization = FALSE, columnWise = TRUE, biomks = svars,  
                    outcomeType = outcomeType, Y = Y, time = time, event = event, outfile = paste0("LASSO_reg_",outfile))
  }
  return(fit)
}
