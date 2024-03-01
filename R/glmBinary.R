#' A function rather aimed at developers
#' @importFrom stats glm
#' @importFrom stats as.formula
#' @importFrom stats coef
#' @importFrom stats confint.default
#' @noRd


glmBinary = function(datain, Yin, Xin){
  ## datain is a data frame
  ## Yin is a binary outcome name
  ## Xin is one x variable, or a group of x variables
  f1 =  as.formula(paste(Yin, paste(Xin, collapse=" + "), sep=" ~ "))
  
  fit1 = glm(formula = f1, family="binomial", data = datain)  
  
  out1 = summary(fit1)$coef
  out2 = exp(cbind("Odds_ratio" = coef(fit1), confint.default(fit1, level = 0.95)))
  
  out12 = cbind(out1,out2)
  return(list(fit1,out12))
}

