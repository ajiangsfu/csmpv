#' rms_model is to use rms package to build up a predictive model
#' @description This is the function to build up a predictive model based on rms with variables selected from LASSO_plus or a model from other R programs/fucntions.
#' @details This function is to build up a predictive model after variable selection with LASSO_plus or a model from other R programs/fucntions. 
#'       R package rms allows us to easily generate calibration (bootstrap) and nomograph plots, calculate C-index when necessary, predict outcome for a new data set. 
#'       If a rms based model is not needed, the model from LASSO_plus object can be directly used as well.
#'       For both of these models, a new data set should be comparable to the data set used for variable selection
#'       rms::calibrate, rms::nomogram, rms::rcorr.cens, rms::Predict
#' @param afit a LASSO_plus object returned by LASSO_plus, or any other R functions/programs
#' @param data a data frame that used to get afit, only need it if the outcome variable is time to event, otherwise data is extracted from afit
#' @param newdata a new data frame for prediction
#' @param outfile a string to indicate output file name without file type but should include all path information
#' @keywords predictive model
#' @author Aixiang Jiang
#' @return A list with three items is returned: 
#' \item{fit}{a model with selected variables for the given outcome variable}
#' \item{coefs}{model coefficients and 95% CI}
#' @references 
#' 
#' @export
#' 
#' wrote code first, later, update the above section

rms_model = function(afit, data = NULL, newdata = NULL, outfile){

  if(is.null(afit)){
    stop("Please input a model object")
  }
  #data = afit$data# realize that coxph does not output data
  # I could return data from LASSO_plus, then I need to modify the parameters if a regular fit object is in instead of LASSO_plus object
  # if I do not want to change on LASSO_plus, then I need to ask for data input if the model type is cox
  
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
  
  ## I should use a list to return all I need, and save 
  newfit = NULL
  
  if(coxtype){
    # leave the common steps outside of if - statements, only work on specific steps I need for a specific outcome
    tmp = strsplit(toString(aform), split = "Surv\\(")
    tmp = strsplit(tmp[[1]][2], split = ",")[[1]][1]
    units(data[,tmp]) = "Year"  ## will set "Year" as parameter later
    newfit = rms::cph(formula = aform, data = data, x=TRUE, y=TRUE, surv = T)
    
    # use parse to convert a string to an expression
    # tmp = strsplit(toString(aform), split = ",")[[1]][2:3] ## this is tricky, the order is changed after toString
    # tmp = toString(tmp)
    # surv.obj=with(data,tmp) #This will be used for rcorr.cens

    ###Create your survival estimates
    estimates=rms::survest(newfit,newdata=data,times=2)$surv  ## use 2 years for now, will make it as a parameter later
    
    ###Determine concordance
    # vcindex = Hmisc::rcorr.cens(x=estimates,S=surv.obj)
    ## it is too bad, the above surv.obj still does not work out
    
    ## try again:
    vcindex = Hmisc::rcorr.cens(x=estimates,S=newfit$y) 
    ## works now, this is a vector!
    
    write.csv(vcindex, paste0(outfile, "Cindex.csv"))
    ## the calibration plot and nomograph should be put here
    
    calib = rms::calibrate(newfit, B= 200, u = 2) ## again, use 2 years as example, will treat it as a parameter later
    
    pdf(paste0(outfile, "_bootstrap_calibration.pdf"))
    plot(calib)
    dev.off()
 
    #survobj = rms::Survival(newfit)
    # pdf(paste0(outfile, "_summary_nomograph.pdf"))
    # plot(summary(newfit), log = T)
    # plot(rms::nomogram(newfit,fun=function(x) survobj(2, x), funlabel="2-year"), cex.axis = 0.7)
    # ## again, use two years as example, but should set it as a parameter
    # dev.off()
    ## maybe it is not necessary to estimate 2 years survival prob, thus, do it all together in the end

  }else if(glmtype){
    newfit = rms::lrm(formula = aform, data = data, x=TRUE, y=TRUE)
   #  stats	
   #  vector with the following elements: number of observations used in the fit, maximum absolute value of first derivative of 
   # log likelihood, model likelihood ratio chi-square, d.f., P-value, c index (area under ROC curve), Somers' D_{xy}, 
   # Goodman-Kruskal gamma, Kendall's tau-a rank correlations between predicted probabilities and observed response, 
   # the Nagelkerke R^2 index, the Brier score computed with respect to Y > its lowest level, the g-index, 
   # gr (the g-index on the odds ratio scale), and gp (the g-index on the probability scale using the same cutoff used 
   # for the Brier score). Probabilities are rounded to the nearest 0.0002 in the computations or rank correlation indexes. 
   # In the case of penalized estimation, the "Model L.R." is computed without the penalty factor, and "d.f." is the effective d.f. 
   # from Gray's (1992) Equation 2.9. The P-value uses this corrected model L.R. chi-square and corrected d.f. 
   # The score chi-square statistic uses first derivatives which contain penalty components.
    # Obs    Max Deriv   Model L.R.         d.f.            P            C          Dxy        Gamma        Tau-a           R2 
    # 4.900000e+01 2.303972e-05 3.577563e+01 6.000000e+00 3.047679e-06 9.431373e-01 8.862745e-01 8.862745e-01 3.843537e-01 7.315664e-01 
    # R2(49)     R2(6,49)     R2(31.2)   R2(6,31.2)        Brier            g           gr           gp 
    # 5.181464e-01 4.553795e-01 6.820164e-01 6.146482e-01 7.827159e-02 7.231361e+00 1.382102e+03 3.867713e-01 
    ## C index is included, save all
    
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
  
  #calib = rms::calibrate(newfit, B= 200, u = 2) ## again, use 2 years as example, will treat it as a parameter later
  # for time to event, I do need to set u: time point, and set units of it, so, cannot do it here
  # pdf(paste0(outfile, "_bootstrap_calibration.pdf"))
  # plot(calib)
  # dev.off()
  
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


############# testing code, will remove later ###########
## read in data first:
library(rstudioapi)
current_working_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
# "/Users/aijiang/Desktop/AJworking/SepOct2022/varSelectPred_SepOct2022/varSelectPred_package/varSelectPred/R"
setwd(current_working_dir)
dat = read.csv("../../exampleData/proteinPerc_ASCT1_postBMTFFS_49ids_20220331.csv", header = T, row.names = 1, stringsAsFactors = F)
tmp = grep("Percent$", colnames(dat))
vars = colnames(dat)[tmp]


## binary
bfit = LASSO_plus(data = dat, biomks = vars, Y =  "CODE_postBMTFFS", outfile = "binaryOutcome_20221114")

## continuous
cfit = LASSO_plus(data = dat, biomks = vars, outcomeType = "continuous", Y =  "postBMTFFS", outfile = "continousOutcome_20221114")

## time to event
tfit = LASSO_plus(data = dat, biomks = vars, outcomeType = "time-to-event", time = "postBMTFFS",  event = "CODE_postBMTFFS", 
           outfile = "timeToEvenOutcome_20221114")

## all I need: bfit$data, and bfit$formula, and bfit$call
## write code here, later put it into the rms_model function
# > bfit$call
# glm(formula = as.formula(paste(Y, vX, sep = " ~ ")), family = "binomial", 
#     data = data)
# > cfit$call
# glm(formula = as.formula(paste(Y, vX, sep = " ~ ")), data = data)
# > tfit$call
# survival::coxph(formula = as.formula(paste(survY, survX, sep = " ~ ")), 
#                 data = data)

## everything works
## I should:
## 1) remove the results file to my presentation folder together with this rms_model.R
## 2) post to GitHub
## 3) clean up all current files, build up package, and re-post
## 4) work on PPT file
## 5) if I still have time, I should work on the last part: validation







