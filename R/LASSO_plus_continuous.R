#' A function rather aimed at developers
#' @noRd
#' 
#' 
#' 
#' write it as an internal function for now, may change it later to be a function that can be used with help file

LASSO_plus_continuous = function(data, biomks,  Y, topN = 10){
  
  vars = intersect(biomks, colnames(data))

  # ### remove variables if their sd < 0.0000001
  sdcol = apply(data[, vars],2,sd, na.rm = TRUE)
  sdtmp = which(sdcol>0.0000001)
  sdtmp = names(sdcol)[sdtmp]
  data = data[,c(Y, sdtmp)] 
  
  ### remove NA for Y side if there are any
  natmp = which(is.na(data[, Y]))
  
  if(length(natmp) > 0){
    data = data[-natmp,]
  }

  lassoF = glmnet::glmnet(x= data.matrix(data[,vars]), y=as.numeric(data[,Y]))
  
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
  
  ### now, let's try to use single variable approach
  #survY = paste0("Surv(", morevar[1],",", morevar[2], ")")  ## this is duplicated, and should remove
  bout = t(sapply(vars, function(var1){
    out1 = glm(as.formula(paste(Y, var1, sep=" ~ ")), data = data)
    return(summary(out1)$coefficients[2,])
  }))
  
  colnames(bout) = c("coef", "se(coef)", "z", "p")
  ### get FDR
  bout = data.frame(bout)
  bout$adjusteP = p.adjust(bout$p, method = "BH")
  
  ### order them and write out
  bout = bout[order(bout$adjusteP),]
  ## cutoff, FDR < 0.05
  bout = subset(bout, bout$adjusteP <= 0.05)
  # if(dim(bout)[1] > 0){ 
  #   bout = head(bout, n = topN)
  #   write.csv(bout, gsub("LASSO", "singleVar", lassoCoef))
  # }  
  topn2 = rownames(bout)  ### even this is empty, the following union still works
  
  ### combine topn1 and topn2, run a full model and step
  topn1 = union(topnsh,topn2)
  
  if(length(topn1) < 1){
    stop("No significant varaible was selected by LASSO_plus")
  }
  
  vX = paste(topn1, collapse=" + ")
  
  fit1 = glm(as.formula(paste(Y, vX, sep=" ~ ")), data = data)
  
  ###########################################################
  ####### the following are LASSO_plus outputs ################
  
  fit2 = step(fit1)
  
  fitout = summary(fit2)$coef
  
  # ### remove intercept
  fitout = fitout[-1,]
  
  ### re-run fit1 based on variables in fit2
  vX = paste(rownames(fitout), collapse=" + ")
  fit =  glm(as.formula(paste(Y, vX, sep=" ~ ")), data = data)
  
  coefs = summary(fit1)$coef
  
  return(list(fit,coefs))
  
}


### I put the R code file and the data file under the same folder:
library(rstudioapi)
current_working_dir = dirname(rstudioapi::getActiveDocumentContext()$path)

setwd(current_working_dir)

dat = read.csv("../../exampleData/proteinPerc_ASCT1_postBMTFFS_49ids_20220331.csv")

tmp = grep("_Percent$", colnames(dat))
vars = colnames(dat)[tmp]

out = LASSO_plus_continuous(data = dat,  biomks = vars,  Y = "postBMTFFS", topN = 10)

