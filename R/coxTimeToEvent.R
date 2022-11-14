#' A function rather aimed at developers
#' @noRd

coxTimeToEvent = function(datain, time, event, Xin){
  ## datain is a data frame
  ## the outcome type is time to event, so, there are two outcome variable names: time and event
  ## Xin is one x variable, or a group of x variables

  survY = paste0("survival::Surv(", time,",", event, ")")
  survX = paste(Xin, collapse=" + ")
  
  fit1 = survival::coxph(as.formula(paste(survY, survX, sep=" ~ ")), data = data)
  
  ## the following part is from my old function: coxphOutput
  
  coxphObj = summary(fit1)
  
  ### focus on $coef, $conf.int
  coxcoef = coxphObj$coefficients
  coxci = coxphObj$conf.int
  
  coxout = coxcoef
  
  mm = dim(coxout)[2]
  
  if(dim(coxout)[1] == 1){
    tmp1 = as.numeric(coxout)
    tmp2 = as.numeric(coxci)
    coxout = c(tmp1[1:2], tmp2[3:4], tmp1[(mm-2):mm])
    coxout = data.frame(matrix(coxout, nrow = 1, ncol = 7))
    
  }else{
    coxout = cbind(coxout[,c(1:2)], coxci[,3:4],coxout[,(mm-2):mm])
  }
  colnames(coxout) = c("beta", 'HR', "HR_95%CI_lower","HR_95%CI_upper", "beta_se", "beta_z", "Pvalue_beta")

  return(list(fit1,coxout))
}

