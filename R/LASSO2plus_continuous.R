
#' A function rather aimed at developers
#' @importFrom stats glm
#' @importFrom stats as.formula
#' @importFrom stats p.adjust
#' @importFrom stats step
#' 
#' @noRd

LASSO2plus_continuous = function(data, biomks,  Y, outfile = "nameWithPath"){

  vars = intersect(biomks, colnames(data))
  
  # # ### remove variables if their sd < 0.0000001
  # sdcol = apply(data[, vars],2,sd, na.rm = TRUE)
  # sdtmp = which(sdcol>0.0000001)
  # sdtmp = names(sdcol)[sdtmp]
  # data = data[,c(Y, sdtmp)] 
  
  ### remove NA for Y side if there are any
  natmp = which(is.na(data[, Y]))
  
  if(length(natmp) > 0){
    data = data[-natmp,]
  }
  
  ## only need to get selected variable names from LASSO2 and save it as topnsh
  lasres = LASSO2(data = data, biomks = biomks, outcomeType = "continuous", Y = Y, outfile = outfile)
  topnsh = names(lasres$coefs)
  
  ###############################################################################
  ###################### the above are LASSO2 results ################
  ###############################################################################
  
  ### now, let's try to use single variable approach
  #survY = paste0("Surv(", morevar[1],",", morevar[2], ")")  ## this is duplicated, and should remove
  bout = t(sapply(vars, function(var1){
    out1 = glm(as.formula(paste(Y, var1, sep=" ~ ")), data = data)
    #return(summary(out1)$coefficients[2,])  
    ## this is not quite right, when there are multiple levels for a single categorical, only the 2nd instead of the last level vs the 1st level is used
    ## in order to keep one sig contrast, do the following
    coefout = summary(out1)$coefficients
    if(dim(coefout)[1] == 2){
      coefout = coefout[2,]
    }else{
      tmp = coefout[-1,4]
      tmp = which.min(tmp)
      coefout = coefout[tmp+1,]
    }
    return(coefout)
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
  lastvars = colnames(fit2$model)[-1]
  vX = paste(lastvars, collapse=" + ")
  fit =  glm(as.formula(paste(Y, vX, sep=" ~ ")), data = data)
  
  coefs = summary(fit1)$coef
  
  return(list(fit,coefs))
  
}