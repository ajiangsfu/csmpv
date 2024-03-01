#' A function rather aimed at developers
#' @importFrom stats glm
#' @importFrom stats as.formula

#' @noRd


lmContinuous = function(datain, Yin, Xin){
  ## datain is a data frame
  ## Yin is a continuous outcome name
  ## Xin is one x variable, or a group of x variables
  f1 =  as.formula(paste(Yin, paste(Xin, collapse=" + "), sep=" ~ "))
  
  ## change on 20230718, use glm instead of lm
  fit1 = glm(formula = f1, data = datain)  
  
  out1 = summary(fit1)$coef

  return(list(fit1,out1))
}

