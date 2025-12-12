
#' A function rather aimed at developers
#' @importFrom stats as.formula
#' @importFrom stats p.adjust
#' @importFrom stats as.formula
#' @importFrom stats step
#' 
#' @noRd

LASSO_plus_timeToEvent = function(data, biomks,  time, event, topN = 10){
  
  vars = intersect(biomks, colnames(data))
 
  ### remove NA for outcome side if there are any
  natmp = which(is.na(data[, event]))
  
  if(length(natmp) > 0){
    data = data[-natmp,]
  }
  
  surObj = survival::Surv(data[,time], data[,event])
  
  lassoF = glmnet::glmnet(x= data.matrix(data[,vars]), y= surObj, family = "cox", type.measure = "C")
  allcoefs = data.matrix(lassoF$beta)
  
  ## calculate how many variables are remained, and a0 is in a different output 
  coef01 = apply(allcoefs, c(1,2), FUN = function(xx) { ifelse(abs(xx) > 0.00001, 1, 0)})
  csum = colSums(coef01)
  
  tt = csum[duplicated(csum)]
  
  ### to find who close to top N
  tmp = which.min(abs(tt - topN))  ### it gives the 1st matching
  
  ### now, find lambda for this tmp
  scoefs = allcoefs[,names(tmp)]
  
  ##############################
  ### get out the variable names whose coefs are not 0
  topnsh = which(abs(scoefs) > 0.00001)
  
  ### get variable names only
  topnsh = names(topnsh)  ### this is the variable list from LASSO
  
  if(length(topnsh) < 1){
    stop("No significant varaible was selected by LASSO")
  }
  
  ###############################################################################
  ###################### the above are customized LASSO results ################
  ###############################################################################
  survY = paste0("survival::Surv(", time,",", event, ")")
  cox1out = t(sapply(vars, function(var1){
    out1 = survival::coxph(as.formula(paste(survY, var1, sep=" ~ ")), data = data)
    coefout = summary(out1)$coefficients
    if(dim(coefout)[1] > 1){
      tmp = coefout[,5]
      tmp = which.min(tmp)
      coefout = coefout[tmp,]
    }
    return(coefout)
  }))
  
  colnames(cox1out) = c("coef", "exp(coef)","se(coef)", "z", "p")
  ### get FDR
  cox1out = data.frame(cox1out)
  cox1out$adjusteP = p.adjust(cox1out$p, method = "BH")
  
  ### order them 
  cox1out = cox1out[order(cox1out$adjusteP),]
  ## cutoff, FDR < 0.05
  cox1out = subset(cox1out, cox1out$adjusteP <= 0.05)
  
  topn2 = rownames(cox1out)
  
  ### combine topn1 and topn2, run a full model and step
  topn1 = union(topnsh,topn2)
  
  if(length(topn1) < 1){
    stop("No significant varaible was selected by LASSO_plus")
  }
  
  ### now, work on stepwise 
  vX = paste(topn1, collapse=" + ")
  
  fit1 = survival::coxph(as.formula(paste(survY, vX, sep=" ~ ")), data = data)
  
  fit2 = step(fit1)
 
  ### final model
  lastvars = names(fit2$assign)
  survX = paste(lastvars, collapse=" + ")
  
  # if less than 2 variables in lastvars, use the list before stepwise selection
  if(length(lastvars) < 2){
    lastvars = topn1
  }
  
  fit = survival::coxph(as.formula(paste(survY, survX, sep=" ~ ")), data = data)
  
  coxphObj = summary(fit)
  
  ### focus on $coef, $conf.int
  coxcoef = coxphObj$coefficients
  coxci = coxphObj$conf.int
  
  coefs = coxcoef
  
  mm = dim(coefs)[2]
  
  if(dim(coefs)[1] == 1){
    tmp1 = as.numeric(coefs)
    tmp2 = as.numeric(coxci)
    coefs = c(tmp1[1:2], tmp2[3:4], tmp1[(mm-2):mm])
    coefs = data.frame(matrix(coefs, nrow = 1, ncol = 7))
    
  }else{
    coefs = cbind(coefs[,c(1:2)], coxci[,3:4],coefs[,(mm-2):mm])
  }
  colnames(coefs) = c("beta", 'HR', "HR_95%CI_lower","HR_95%CI_upper", "beta_se", "beta_z", "Pvalue_beta")
  
  return(list(fit,coefs))
}