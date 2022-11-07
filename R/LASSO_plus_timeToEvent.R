#' A function rather aimed at developers
#' @noRd
#' 
#' 
#' 
#' write it as an internal function for now, may change it later to be a function that can be used with help file

LASSO_plus_timeToEvent = function(data, biomks,  time, event, topN = 10){
  
  vars = intersect(biomks, colnames(data))
  
  # ### remove variables if their sd < 0.0000001
  sdcol = apply(data[, vars],2,sd, na.rm = TRUE)
  sdtmp = which(sdcol>0.0000001)
  sdtmp = names(sdcol)[sdtmp]
  data = data[,c(time, event, sdtmp)] 
  
  ### remove NA for outcome side if there are any
  natmp = which(is.na(data[, event]))
  
  if(length(natmp) > 0){
    data = data[-natmp,]
  }

  surObj = survival::Surv(data[,time], data[,event])
  
  lassoF = glmnet::glmnet(x= data.matrix(data[,vars]), y= surObj, family = "cox")
  allcoefs = data.matrix(lassoF$beta)

#  is it necessary to save the lassoplot and lambda tables?
#   outfile = paste0("glm_",topN,"_",filename, "_", title)
#   lassoPlot = paste0(outfile, "_LASSO.pdf")
#   lassoCoef = gsub("pdf", "csv",lassoPlot)
#   
#   pdf(lassoPlot)
#   plot(lassoF)
#   #plotres(lassoF)
#   plot(lassoF,xvar="lambda",label=TRUE)
#   dev.off()
  
  ####################################################

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
    return(summary(out1)$coefficients)
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
  afit1 = summary(fit2)$coef
  survX = paste(rownames(afit1), collapse=" + ")
  
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



### I put the R code file and the data file under the same folder:
library(rstudioapi)
current_working_dir = dirname(rstudioapi::getActiveDocumentContext()$path)

setwd(current_working_dir)

dat = read.csv("../../exampleData/proteinPerc_ASCT1_postBMTFFS_49ids_20220331.csv")

tmp = grep("_Percent$", colnames(dat))
vars = colnames(dat)[tmp]

out = LASSO_plus_timeToEvent(data = dat,  biomks = vars, time = "postBMTFFS",  event = "CODE_postBMTFFS", topN = 10)

