
#' A Wrapper for Building Predictive Models using the rms Package
#' 
#' @description This wrapper function harnesses the capabilities of the rms package to construct robust predictive models using various modeling techniques.
#' Whether you're working with linear regression (lm), generalized linear models (glm), or Cox proportional hazards models (coxph) directly or generating these 
#' objects from other functions, this function streamlines the model-building process and enhances analysis. The function features built-in bootstrap-based 
#' internal validation, model score generation, and result storage.
#'
#' @details The wrapper function capitalizes on the power of the rms package, providing a seamless interface for creating predictive models.
#' With the rms package, the function conducts efficient bootstrap-based internal validation, enabling thorough model evaluation. 
#' Nomograph plots and C-index calculations are also at your disposal. Furthermore, the function adapts to your needs by predicting model scores for new datasets. 
#' In the absence of new data, it generates predictions based on the training data. However, when you provide an independent dataset complete with outcome 
#' variable information, the function extends to external validation. This enables assessing model performance in an independent context using the rms R package.
#'
#' All results, from the model itself to internal validation metrics and predictions, are conveniently stored locally for easy access and further analysis.

#' @import rms
#' @importFrom stats anova
#' @importFrom utils write.csv
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' 
#' @param afit A model fit from lm, glm, or coxph.
#' @param data A data frame used to obtain the 'afit' object. This parameter is only required when the outcome variable is time-to-event; otherwise, the data is extracted from 'afit'.
#' @param newdata A new data frame for prediction.
#' @param newY A logical variable indicating whether 'newdata' contains the outcome variable.
#' @param u A single numeric value representing follow-up time for survival outcomes, used to estimate the survival probability for the given time point.
#' @param outfile A string indicating the output file name without the file type extension, but including the complete path information.
#' @return Predictions are returned.

#' @author Aixiang Jiang
#' @references 
#'   Harrell Jr F (2023). rms: Regression Modeling Strategies. R package version 6.7-1, <https://CRAN.R-project.org/package=rms>
#'   
#'   Harrell Jr F (2023). Hmisc: Harrell Miscellaneous. R package version 5.1-1, <https://CRAN.R-project.org/package=Hmisc>

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
#' # Here, we use continuous as an example:
#' clr = LASSO2_reg(data = tdat, biomks = Xvars,
#'                  outcomeType = "continuous", Y = "Age",
#'                  outfile = paste0(temp_dir, "/continuousLASSO2_reg"))
#' pclr = rms_model(clr$fit, newdata = vdat,
#'                  outfile = paste0(temp_dir, "/pred_LASSO2reg_continuous"))
#' # You might save the files to the directory you want.
#' 
#' # To delete the "temp_dir", use the following:
#' unlink(temp_dir)

#' @export

rms_model = function(afit, data = NULL, newdata = NULL, newY = FALSE, u = 2, outfile){
  
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
  
  if(!glmtype & !coxtype){
    lmtype = TRUE
  }
  
  if(is.null(data) & !coxtype){
    data = afit$data
  }else if (is.null(data) & coxtype){
    stop("Please input a data set")
  }
   
  dd = rms::datadist(data)
  options(datadist= dd)

  newfit = NULL
  
  ## maybe it is a good idea to save objects that I need for plots
  bootvalidate = NULL # this is for training data set
  newplot = NULL # this is for new data set

  if(coxtype){
    # leave the common steps outside of if - statements, only work on specific steps I need for a specific outcome
    tmp = strsplit(toString(aform), split = "Surv\\(")
    tmpt = strsplit(tmp[[1]][2], split = ",")[[1]][1]
    tmpe = strsplit(tmp[[1]][2], split = ",")[[1]][2]
    tmpe = gsub("\\)", "", tmpe)
    tmpe = gsub(" ", "", tmpe)
    
    units(data[,tmpt]) = "Year"  ## will set "Year" as parameter later
    newfit = rms::cph(formula = aform, data = data, x=TRUE, y=TRUE, surv = T)
    
    estimates=rms::survest(newfit,newdata=data,times=u)$surv
    ## this is survival prob for all in newdata at a given time point: u
    
    cindex = Hmisc::rcorr.cens(x=estimates,S=newfit$y) 
    write.csv(cindex, paste0(outfile, "training_Cindex.csv"))
    
    bootvalidate = rms::validate(newfit, B= 200)
    Cindex = bootvalidate[1,1:5]/2 + 0.5
    bootvalidate = rbind(bootvalidate, c(Cindex, bootvalidate[1,6]))
    nv = dim(bootvalidate)[1]
    rownames(bootvalidate)[nv] = "Cindex"
    ## for the optimism, I need -0.5
    bootvalidate[nv, 4] = bootvalidate[nv, 4] - 0.5
 
    ## independent data validation when newdata contains outcome variable
    if(!is.null(newdata)){
      # note on 20230712, To estimate survival probability, I do need to define a follow-up year; however, 
      # however, if newdata contain time variable, I should use this variable instead of a constant follow-up time point 
      #     for survival probability estimation
      
      ## without newY
      ss = survest(newfit,newdata=newdata,times=u)$surv
      survs = Surv(time = newdata[,tmpt], event = newdata[,tmpe])
      cindex = rcorr.cens(x=ss, S = survs)
      write.csv(cindex, paste0(outfile, "_newData_Cindex.csv"))
      write.csv(ss, paste0(outfile,"SurvProb_at_year_",u ,"_newData.csv"))
      
      if(newY){
        vout = val.surv(newfit,newdata, S = survs, est.surv = ss)
        
        pdf(paste0(outfile, "_newData_validation.pdf"))
        plot(vout, xlab = "Predicted Probability" , 
                      ylab = "Actual Probability", main = paste0("C-index: ", format(cindex[1], digits = 4)))
        dev.off()
        # 
        plot(vout, xlab = "Predicted Probability" , 
             ylab = "Actual Probability", main = paste0("C-index: ", format(cindex[1], digits = 4)))

      }
      
    }

  }else if(glmtype){
    newfit = rms::lrm(formula = aform,data = data, x=TRUE, y=TRUE)
    write.csv(newfit$stats, paste0(outfile, "training_Cindex.csv"))
    
    bootvalidate = rms::validate(newfit, B= 200) 
    Cindex = bootvalidate[1,1:5]/2 + 0.5
    bootvalidate = rbind(bootvalidate, c(Cindex, bootvalidate[1,6]))
    rownames(bootvalidate)[dim(bootvalidate)[1]] = "Cindex"
    nv = dim(bootvalidate)[1]
    rownames(bootvalidate)[nv] = "Cindex"
    ## for the optimism, I need -0.5
    bootvalidate[nv, 4] = bootvalidate[nv, 4] - 0.5

    if(!is.null(newdata)){
      ## val.prob provide both newdata validation plot and c-index
      tmp = strsplit(toString(aform), split = ",")
      tmp = tmp[[1]][2]
      tmp = gsub(" ", "", tmp) ## finally we get Y var name
      ## the example: for "counts ~ outcome + treatment", toString of it: "~, counts, outcome + treatment"
      pred.logit = stats::predict(newfit, newdata = newdata) ## this is model score
      phat = 1/(1+exp(-pred.logit)) # this is prob
      write.csv(phat, paste0(outfile,"Prob_posGroup_newData.csv"))
      
      if(newY){
        pdf(paste0(outfile, "_newData_validation.pdf"))
        rms::val.prob(phat, newdata[,tmp])
        dev.off()
        
        rms::val.prob(phat, newdata[,tmp])
      }
    }
    
  }else if(lmtype){
    newfit = rms::ols(formula = aform,data = data, x=TRUE, y=TRUE)
    ## no c-index, but there are some info about model including R**2, LR: likelihood ratio ch-square stat
    write.csv(newfit$stats, paste0(outfile, "training_R2_LR.csv"))
    
    bootvalidate = rms::validate(newfit, B= 200) 
 
    ## for validation: maybe simply use a regression line with RMSE?
    if(!is.null(newdata)){
      tmp = strsplit(toString(aform), split = ",")
      tmp = tmp[[1]][2]
      tmp = gsub(" ", "", tmp) ## finally we get Y var name
      yhat = stats::predict(newfit, newdata = newdata) ## this is model score
      write.csv(yhat, paste0(outfile,"predicted_modelScore_newData.csv"))
      
      if(newY){
        newplot = validate_continuous(yhat = yhat, yobs = newdata[,tmp])
        pdf(paste0(outfile, "_newData_validation.pdf"))
        print(newplot)
        dev.off()
        print(newplot)
      }
    }   
  }
  
  if(is.null(newdata)){
    newdata = data
  }
  prediction = stats::predict(newfit, newdata = newdata) ## when new data is NULL, use the current data for "prediction"

  modelout = data.frame(anova(newfit))
  validat = rms::validate(newfit, B = 200)
  
  logflag = TRUE
  if (lmtype){
    logflag = FALSE
  }
  
  pdf(paste0(outfile, "_trainingPlots_.pdf"))
  plot(summary(newfit), log = logflag)
  plot(rms::nomogram(newfit), cex.axis = 0.6)
  dev.off()
  
  plot(summary(newfit), log = logflag)
  plot(rms::nomogram(newfit), cex.axis = 0.6)
  print(bootvalidate)
 
  ## write out: prediction, modelout, and validation
  write.csv(prediction, paste0(outfile,"_modelScore.csv"))
  write.csv(modelout, paste0(outfile,"_rms_model.csv"))
  write.csv(validat, paste0(outfile,"_bootstrap_validation.csv"))
  
  return(prediction)

}
