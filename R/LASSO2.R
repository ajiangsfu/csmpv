#' Variable Selection using Modified LASSO with a Minimum of Two Remaining Variables
#' @description
#' This function conducts variable selection using LASSO (Least Absolute Shrinkage and Selection Operator) with a minor adaptation. 
#' It calculates the mean lambda value from multiple cv.glmnet runs and ensures the selection of at least two variables.
#'
#' @details
#' The function utilizes glmnet::cv.glmnet for cross-validation-based variable selection with 
#' the largest value of lambda such that error is within 1 standard error of the minimum. 
#' To mitigate randomness from cross-validation splits, it conducts 10 runs (this number can later be parameterized) of n-fold cv.glmnet. The resulting 
#' average lambda value across these runs serves as the final lambda. Subsequently, the final regularization regression is performed on the 
#' complete dataset using this mean lambda value. Following this, the function assesses the count of remaining variables. If only one or none 
#' are selected, the function defaults to selecting the first lambda that results in at least two chosen variables on the full dataset.
#' This function is designed to handle three types of outcome variables: continuous, binary, and time-to-event.
#' @param data A data matrix or a data frame where samples are in rows and features/traits are in columns.
#' 
#' @param standardization A logical variable indicating if standardization is needed before variable selection. The default is FALSE.
#' @param columnWise A logical variable indicating if column-wise or row-wise normalization is needed. The default is TRUE, which means column-wise normalization is performed.
#' This is only meaningful when "standardization" is TRUE.
#' @param biomks A vector of potential biomarkers for variable selection. They should be a subset of the column names in the "data" variable.
#' @param outcomeType The outcome variable type. There are three choices: "binary" (default), "continuous", and "time-to-event".
#' @param Y The outcome variable name when the outcome type is either "binary" or "continuous".
#' @param time The time variable name when the outcome type is "time-to-event".
#' @param event The event variable name when the outcome type is "time-to-event".
#' @param nfolds The number of folds for cross-validation. The default is 10.
#' @param outfile A string representing the output file, including the path if necessary, but without the file type extension.
#' @import glmnet
#' @importFrom Matrix nnzero
#' @importFrom grDevices pdf
#' @importFrom graphics title
#' @importFrom stats coef
#' @importFrom grDevices dev.off
#' 
#' @author Aixiang Jiang
#' @return A list is returned:
#' \item{coefs}{A vector of LASSO coefficients}
#' \item{h0}{Cumulative baseline hazard table, for time to event outcome only}
#' \item{Y}{The outcome variable name when the outcome type is either "binary" or "continuous".}
#' \item{time}{The time variable name when the outcome type is "time-to-event".}
#' \item{event}{The event variable name when the outcome type is "time-to-event".}
#' \item{standardization}{A logical variable indicating if standardization is needed before variable selection.}
#' \item{columnWise}{A logical variable indicating if column-wise or row-wise normalization is needed.}
#' \item{outcomeType}{The outcome variable type.}
#' \item{allplot}{A plot object}
#' A shrunken coefficient vector is returned 
#' @references 
#'   Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent (2010), Journal of Statistical Software, Vol. 33(1), 1-22, doi:10.18637/jss.v033.i01.
#'   
#'   Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol. 39(5), 1-13, doi:10.18637/jss.v039.i05.
#'
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
#' # Here, we use time-to-event as an example:
#' # tl = LASSO2(data = tdat, biomks = Xvars,
#' #             outcomeType = "time-to-event",
#' #             time = "FFP..Years.",event = "Code.FFP",
#' #             outfile = paste0(temp_dir, "/survivalLASSO2"))
#' # You might save the files to the directory you want.
#' 
#' # To delete the "temp_dir", use the following:
#' unlink(temp_dir)


#' @export

LASSO2 = function(data = NULL, standardization = FALSE, columnWise = TRUE, biomks = NULL, outcomeType = c("binary","continuous","time-to-event"), 
                      Y = NULL, time = NULL, event = NULL, nfolds = 10, outfile = "nameWithPath"){
  
  if(is.null(data)){
    stop("Please input a data set")
  }
  if(standardization){
    data = standardize(data, byrow = !columnWise)
  }
  
  alls = NA
  h0 = NA
  outcomeType = outcomeType[1]
  
  vars = intersect(biomks, colnames(data))
  
  ## note on 20231123 about type.measure
  ## 1) for binary, use default type.measure="deviance"
  ## 2) for continuous, use default type.measure="deviance"
  ## 3) for cox, use type.measure="C"
  
  
  ## add more on 20231024: when there is less than 2 selected variables, set it to at least 2
  if(outcomeType == "binary"){
    
    ## change to multiple runs on 20231122 to reduce randomness caused by random splits in cv.glmnet
    lvalue = rep(0,10)
    for (i in 1:10){
      alls = glmnet::cv.glmnet(x= data.matrix(data[,vars]), y=as.numeric(data[,Y]), family = "binomial", nfolds = nfolds)
      # analyze examples in cv.glmnet's help file, lambda.1se is actually used.
      lvalue[i] = alls$lambda.1se
    }
    pdf(paste0(outfile,"_LASSO_output.pdf"))
    plot(alls) 
    title(paste0("Outcome type: ",outcomeType), line = 2.5)
    dev.off()

    # to make sure that model output is always in glmnet object format
    # more change on 20231122, use mean lambda values from mutliple cv.glmnet runs, re-run glmnet with mini lambda after cv.glmnet
    alls = glmnet::glmnet(data.matrix(data[,vars]), y=as.numeric(data[,Y]), family = "binomial", lambda = mean(lvalue))
    nk = nnzero(coef(alls)) - 1 ## for this model, we have intercept, so need to minus 1 to get pure number of selected variables
    
    if(nk < 2){
      tmpres = glmnet::glmnet(data.matrix(data[,vars]), y=as.numeric(data[,Y]), family = "binomial")
      moreres = cbind(tmpres$lambda, tmpres$df)
      min_nonzero_vars <- 2
      selected_lambda <- NULL
      for (i in 1:dim(moreres)[1]){
        nonzero_coefs = moreres[i,2]
        if (nonzero_coefs >= min_nonzero_vars) {
          selected_lambda = moreres[i,1]
          break
        }
      }
      alls = glmnet::glmnet(data.matrix(data[,vars]), y=as.numeric(data[,Y]), family = "binomial", lambda = selected_lambda)
      ## cv.glmenet$glmnet.fit as the same output as glmnet::glmnet
      ## however, in the above, if a specific lambda is given for glmnet::glmnet, it will give a specific output
      ##  to get the same output for a specific lambda for cv.glmenet$glmnet.fit, we need to set the lambda as well
    }
    
  }else if(outcomeType == "continuous"){
    lvalue = rep(0,10)
    for (i in 1:10){
      alls = glmnet::cv.glmnet(x= data.matrix(data[,vars]), y=as.numeric(data[,Y]), nfolds = nfolds)
      # analyze examples in cv.glmnet's help file, lambda.1se is actually used.
      lvalue[i] = alls$lambda.1se
    }
    pdf(paste0(outfile,"_LASSO_output.pdf"))
    plot(alls) 
    title(paste0("Outcome type: ",outcomeType), line = 2.5)
    dev.off()

    # to make sure that model output is always in glmnet object format
    # more change on 20231122, use mean lambda values from mutliple cv.glmnet runs, re-run glmnet with mini lambda after cv.glmnet
    alls = glmnet::glmnet(data.matrix(data[,vars]), y=as.numeric(data[,Y]),lambda = mean(lvalue))
    nk = nnzero(coef(alls)) - 1 ## for this model, we have intercept, so need to minus 1 to get pure number of selected variables
    
    if(nk < 2){
      tmpres = glmnet::glmnet(data.matrix(data[,vars]), y=as.numeric(data[,Y]))
      moreres = cbind(tmpres$lambda, tmpres$df)
      min_nonzero_vars <- 2
      selected_lambda <- NULL
      for (i in 1:dim(moreres)[1]){
        nonzero_coefs = moreres[i,2]
        if (nonzero_coefs >= min_nonzero_vars) {
          selected_lambda = moreres[i,1]
          break
        }
      }
      alls = glmnet::glmnet(data.matrix(data[,vars]), y=as.numeric(data[,Y]), lambda = selected_lambda)
    }
    
  }else if(outcomeType == "time-to-event"){
    surObj = survival::Surv(data[,time], data[,event])
    
    ## change to multiple runs on 20231122 to reduce randomness caused by random splits in cv.glmnet
    lvalue = rep(0,10)
    for (i in 1:10){
      alls = glmnet::cv.glmnet(x= data.matrix(data[,vars]), y= surObj, family = "cox", type.measure = "C", nfolds = nfolds)
      # analyze examples in cv.glmnet's help file, lambda.1se is actually used.
      lvalue[i] = alls$lambda.1se
    }
    pdf(paste0(outfile,"_LASSO_output.pdf"))
    plot(alls) 
    # A plot is produced, and nothing is returned
    # after change to multiple runs, this is actually plotted based on last cv.glmnet, not for all
    #    while this is not much meaningful, at least we get some idea
    title(paste0("Outcome type: ",outcomeType), line = 2.5)
    dev.off()
    # to make sure that model output is always in glmnet object format
    # more change on 20231122, use mean lambda values from mutliple cv.glmnet runs, re-run glmnet with mini lambda after cv.glmnet
    alls = glmnet::glmnet(x= data.matrix(data[,vars]), y= surObj, family = "cox", type.measure = "C",lambda = mean(lvalue))
    nk = nnzero(coef(alls))  ## for this model, no intercept
    
    if(nk < 2){
      tmpres = glmnet::glmnet(x= data.matrix(data[,vars]), y= surObj, family = "cox", type.measure = "C")
      moreres = cbind(tmpres$lambda, tmpres$df)
      min_nonzero_vars <- 2
      selected_lambda <- NULL
      for (i in 1:dim(moreres)[1]){
        nonzero_coefs = moreres[i,2]
        if (nonzero_coefs >= min_nonzero_vars) {
          selected_lambda = moreres[i,1]
          break
        }
      }
      alls = glmnet::glmnet(x= data.matrix(data[,vars]), y= surObj, family = "cox", type.measure = "C",lambda = selected_lambda)
    }
 
    coxfit = coxph(surObj ~ 1)
    h0 = basehaz(coxfit)
  }else{
    stop("Please select the correct outcome type")
  }

  sink(paste0(outfile,"_LASSO_coefs.txt"))
  print(coef(alls))
  sink()
  
  scoefs = coef(alls)
  scoefs = data.frame( predict_names = rownames(scoefs),
                           coef_vals = matrix(scoefs))
  if(scoefs[1,1] == "(Intercept)"){
    scoefs = scoefs[-1,]
  }
  scoefs = subset(scoefs, abs(scoefs[,2]) > 0.0000001)
  coefs = scoefs[,2]
  names(coefs) = scoefs[,1]
  
  outs = list(coefs, h0, Y, time, event, standardization, columnWise, outcomeType, alls) ## when h0=NA, means that this is not time to event outcome, which is OK, same as Y, time, event
  names(outs) = c("coefs", "h0", "Y", "time", "event", "standardization", "columnWise", "outcomeType", "modelObj")
  return(outs) 
}
