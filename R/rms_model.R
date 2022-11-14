#' rms_model is to use rms package to build up a predictive model
#' @description This is the function to build up a predictive model based on rms package after variable selection with LASSO_plus or a model from other R programs/fucntions.
#' @details This function is to build up a predictive model after variable selection with LASSO_plus or a model from other R programs/fucntions. 
#'       R package rms allows us to easily generate calibration (bootstrap) and nomograph plots, calculate C-index when necessary, predict outcome for a new data set. 
#'       If a rms based model is not needed, the model from LASSO_plus object can be directly used as well.
#'       For both of these models, a new data set should be comparable to the data set used for variable selection.
#'       If a new data set is not given, the training data set is used for prediction. In this function, validation and calibration are based on bootstrap.
#' @param afit A LASSO_plus object returned by LASSO_plus, or a model fit from any other R functions/programs
#' @param data A data frame that used to get afit, only this is only needed when the outcome variable is time to event, otherwise data is extracted from afit.
#' @param newdata A new data frame for prediction
#' @param outfile A string to indicate output file name without file type but should include all path information
#' @keywords predictive model
#' @author Aixiang Jiang
#'
#' @references 
#' 
#' @export

rms_model = function(afit, data = NULL, newdata = NULL, outfile){

  if(is.null(afit)){
    stop("Please input a model object")
  }

  aform = afit$formula
  atype = afit$call

  glmtype = grep("binomial", atype)
  coxtype = grep("coxph", atype)
  lmtype = FALSE
  
  if(length(glmtype) > 0){
    glmtype = TRUE
  }else{
    glmtype = FALSE
  }
  
  if(length(coxtype) > 0){
    coxtype = TRUE
  }else{
    coxtype = FALSE
  }
  
  if (!coxtype & !glmtype){
    lmtype = TRUE
  }

  if(is.null(data) & !coxtype){
    data = afit$data
  }else if (is.null(data) & coxtype){
    stop("Please input a data set")
  }
  
  dd = rms::datadist(data)
  options(datadist='dd')
  
  newfit = NULL
  
  if(coxtype){
    # leave the common steps outside of if - statements, only work on specific steps I need for a specific outcome
    tmp = strsplit(toString(aform), split = "Surv\\(")
    tmp = strsplit(tmp[[1]][2], split = ",")[[1]][1]
    units(data[,tmp]) = "Year"  ## will set "Year" as parameter later
    newfit = rms::cph(formula = aform, data = data, x=TRUE, y=TRUE, surv = T)
    
    estimates=rms::survest(newfit,newdata=data,times=2)$surv  ## use 2 years for now, will make it as a parameter later
    vcindex = Hmisc::rcorr.cens(x=estimates,S=newfit$y) 
    write.csv(vcindex, paste0(outfile, "Cindex.csv"))

    calib = rms::calibrate(newfit, B= 200, u = 2) ## again, use 2 years as example, will treat it as a parameter later
    
    pdf(paste0(outfile, "_bootstrap_calibration.pdf"))
    plot(calib)
    dev.off()
 
  }else if(glmtype){
    newfit = rms::lrm(formula = aform, data = data, x=TRUE, y=TRUE)
    write.csv(newfit$stats, paste0(outfile, "Cindex.csv"))
    
    calib = rms::calibrate(newfit, B= 200) ## again, use 2 years as example, will treat it as a parameter later
    
    pdf(paste0(outfile, "_bootstrap_calibration.pdf"))
    plot(calib)
    dev.off()
    
  }else if(lmtype){
    newfit = rms::ols(formula = aform, data = data, x=TRUE, y=TRUE)
    ## no c-index, but there are some info about model including R**2, LR: likelihood ratio ch-square stat
    write.csv(newfit$stats, paste0(outfile, "R2_LR.csv"))
    
    calib = rms::calibrate(newfit, B= 200) ## again, use 2 years as example, will treat it as a parameter later
    
    pdf(paste0(outfile, "_bootstrap_calibration.pdf"))
    plot(calib)
    dev.off()
  }

  if(is.null(newdata)){
    newdata = data
  }
  prediction = predict(newfit, newdata = newdata) ## when new data is NULL, use the current data for "prediction"
  ## this is the predicted model score
  
  modelout = data.frame(anova(newfit))
  validation = rms::validate(newfit, B = 200)
  
  logflag = TRUE
  if (lmtype){
    logflag = FALSE
  }
  
  pdf(paste0(outfile, "_summary_nomograph.pdf"))
  plot(summary(newfit), log = logflag)
  plot(rms::nomogram(newfit), cex.axis = 0.6) ## 0.6 should be set a parameter later
  dev.off()
  
  ## return or write out: prediction, modelout, and validation
  ## maybe just write out:
  write.csv(prediction, paste0(outfile,"_predicted_modelScore.csv"))
  write.csv(modelout, paste0(outfile,"_rms_model.csv"))
  write.csv(validation, paste0(outfile,"_validation_bootstrap.csv"))
}

