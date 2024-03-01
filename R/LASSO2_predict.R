#' Predict and Validate LASSO2 Model Scores
#' 
#' @description
#' This function predicts model scores, positive group probabilities, and risk scores for different outcome types based on a LASSO2 object and a new data set.
#' It also performs validation of predictions when the true outcome variable is available.
#'
#' @param lassoObj A LASSO2 object.
#' @param newdata A new data matrix or data frame where samples are in rows and features/traits are in columns.
#' @param newY The outcome variable for the new data set.
#' @param u A single numeric follow-up time for survival outcomes.
#' @param outfile A string representing the output file, including the path if necessary, but without the file type extension.
#' @author Aixiang Jiang
#' @import glmnet
#' @return A vector of predicted values is returned.
#' @references 
#'   Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent (2010), Journal of Statistical Software, Vol. 33(1), 1-22, doi:10.18637/jss.v033.i01.
#'   
#'   Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol. 39(5), 1-13, doi:10.18637/jss.v039.i05.
#'   
#'   Harrell Jr F (2023). rms: Regression Modeling Strategies_. R package version 6.7-1, <https://CRAN.R-project.org/package=rms>
#'   
#'   Harrell Jr F (2023). Hmisc: Harrell Miscellaneous_. R package version 5.1-1, <https://CRAN.R-project.org/package=Hmisc>
#'
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

#' # The function can work with three different outcome types. 
#' # Here, we use time-to-event as an example:
#' # tl = LASSO2(data = tdat, biomks = Xvars,
#' #             outcomeType = "time-to-event",
#' #             time = "FFP..Years.",event = "Code.FFP",
#' #             outfile = paste0(temp_dir, "/survivalLASSO2"))
#' # To predict the model in a new data set:
#' # ptl = LASSO2_predict(tl, newdata = vdat,
#' #                     outfile = paste0(temp_dir, "/pred_LASSO2_time_to_event"))
#' # You might save the files to the directory you want.
#' 
#' # To delete the "temp_dir", use the following:
#' unlink(temp_dir)

#' @export

LASSO2_predict = function(lassoObj, newdata = NULL,newY = FALSE, u=2, outfile = "nameWithPath"){
  
  coefs = lassoObj$coefs
  
  if(length(coefs) == 0){
    stop("No variable is selected by LASSO")
  }
  baseHz = lassoObj$h0
  Y = lassoObj$Y
  time = lassoObj$time
  event = lassoObj$event
  standardization = lassoObj$standardization
  columnWise = lassoObj$columnWise
  outcomeType = lassoObj$outcomeType

  if(is.null(newdata)){
    stop("Please input a data set")
  }
  if(standardization){
    newdata = standardize(newdata, byrow = !columnWise)
  }

  vars = names(coefs)
  vars = names(coefs)
  
  if(length(vars)>=1 & vars[1] == "(Intercept)"){
    vars = setdiff(vars, "(Intercept)")
  }
  
  model_score = 0 ## since lassobj contains different items for different outcome, models_score will be calculated within if loop

  # ## make some basic transformation and call validation function for each outcome type
  if(outcomeType == "continuous"){
    fullvars = lassoObj$modelObj$beta@Dimnames[[1]]
    model_score = stats::predict(object = lassoObj$modelObj, newx = data.matrix(newdata[,fullvars]))
    model_score = model_score[,1]
    
    ## validation step: only when newY = true, otherwise, return: outs = list(model_score)
    names(model_score) = rownames(newdata)
    outs = model_score
    if(newY){
      yy = newdata[,Y]
      outs = validation(predicted = model_score, outcomeType = "continuous", trueY = yy, outfile = outfile)
    }

  }else if(outcomeType == "binary"){
    fullvars = lassoObj$modelObj$beta@Dimnames[[1]]
    model_score = stats::predict(object = lassoObj$modelObj, newx = data.matrix(newdata[,fullvars]))
    model_score = model_score[,1]
    
    pos_prob = 1/(1+exp(-model_score))
    names(pos_prob) = rownames(newdata)
    outs = pos_prob
    if(newY){
      yy = newdata[,Y]
      outs = validation(predicted = pos_prob, outcomeType = "binary", trueY = yy, outfile = outfile)
    }
  }else if(outcomeType == "time-to-event"){
    fullvars = lassoObj$modelObj$beta@Dimnames[[1]]
    model_score = stats::predict(object = lassoObj$modelObj, newx = data.matrix(newdata[,fullvars]))
    model_score = model_score[,1]
    
    risk_score = exp(model_score)
    names(risk_score) = rownames(newdata)
    outs = risk_score
    if(newY){
      ttime = newdata[, time]
      eevent = newdata[, event]
      outs = validation(predicted = risk_score, outcomeType = "time-to-event", time = ttime, trueEvent = eevent,
                        baseHz = baseHz, u=u, outfile = outfile)
    }
  }
  return(outs)


}





