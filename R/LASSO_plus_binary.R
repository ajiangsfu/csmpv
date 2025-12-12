
#' A function rather aimed at developers
#' @importFrom stats glm
#' @importFrom stats as.formula
#' @importFrom stats p.adjust
#' @importFrom stats step
#' @importFrom stats coef
#' @importFrom stats confint.default
#' @noRd
#' 

LASSO_plus_binary = function(data, biomks,  Y, topN = 10){
  
  vars = intersect(biomks, colnames(data))

  ### remove NA for Y side if there are any
  natmp = which(is.na(data[, Y]))
  
  if(length(natmp) > 0){
    data = data[-natmp,]
  }
  
  lassoF = glmnet::glmnet(x= data.matrix(data[,vars]), y=as.integer(data[,Y]), family = "binomial")
  
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
  
  ### now, let's try to use single variable approach
  bout = t(sapply(vars, function(var1){
    out1 = glm(as.formula(paste(Y, var1, sep=" ~ ")), family="binomial", data = data)
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
  
  topn2 = rownames(bout)  ### even this is empty, the following union still works
  
  ### combine topn1 and topn2, run a full model and step
  topn1 = union(topnsh,topn2)
  
  if(length(topn1) < 1){
    stop("No significant varaible was selected by LASSO_plus")
  }

  vX = paste(topn1, collapse=" + ")
  
  fit1 = glm(as.formula(paste(Y, vX, sep=" ~ ")), family="binomial",data = data)
  
  ###########################################################
  ####### the following are LASSO_plus outputs ################
  
  fit2 = step(fit1)
  
  afit1 = summary(fit2)$coef
  ## I need odds ratio as well
  afit2 = exp(cbind("Odds_ratio" = coef(fit2), confint.default(fit2, level = 0.95)))
  
  tmp = intersect(rownames(afit1), rownames(afit2))
  
  if(length(tmp) > 1){
    fitout = cbind(afit1[tmp,], afit2[tmp,])
  }else{
    stop("No significant varaible was selected by LASSO_plus")
  }
  
  # ### remove intercept, order by odds ratio
  fitout = fitout[-1,]
  
  ### re-run fit1 based on variables in fit2
  lastvars = colnames(fit2$model)[-1]
  
  # if less than 2 variables in lastvars, use the list before stepwise selection
  if(length(lastvars) < 2){
    lastvars = topn1
  }
  
  vX = paste(lastvars, collapse=" + ")
  fit =  glm(as.formula(paste(Y, vX, sep=" ~ ")), family="binomial",data = data)
  
  out1 = summary(fit)$coef
  out2 = exp(cbind("Odds_ratio" = coef(fit), confint.default(fit, level = 0.95)))
  
  coefs = cbind(out1,out2)
  return(list(fit,coefs))

}