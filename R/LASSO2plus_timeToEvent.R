
#' A function rather aimed at developers
#' @importFrom stats as.formula
#' @importFrom stats p.adjust
#' @importFrom stats as.formula
#' @importFrom stats step
#' 
#' @noRd

LASSO2plus_timeToEvent = function(data, biomks,  time, event, outfile = "nameWithPath"){
  
  vars = intersect(biomks, colnames(data))

  ### remove NA for outcome side if there are any
  natmp = which(is.na(data[, event]))
  
  if(length(natmp) > 0){
    data = data[-natmp,]
  }

  ## only need to get selected variable names from LASSO2 and save it as topnsh
  lasres = LASSO2(data = data, biomks = biomks, outcomeType = "time-to-event", time = time, event = event, outfile = outfile)
  topnsh = names(lasres$coefs)
  
  ###############################################################################
  ###################### the above are LASSO2 results ################
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
  
  if(length(names(fit2$coefficients)) < 1){
    stop("No significant varaible was selected by LASSO_plus")
  }
  
  ### final model
  lastvars = names(fit2$assign)
  ## need to make sure the final model contain at least two variables, which is controlled by LASSO2 as well
  ## so, when there is less than 2 variables, just use LASSO2 result
  if(length(lastvars) < 2){
    lastvars = topnsh
  }
  
  survX = paste(lastvars, collapse=" + ")
  
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