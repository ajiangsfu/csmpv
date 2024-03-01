#' A function rather aimed at developers
#' @importFrom stats glm
#' @noRd

gettValue = function(vdat,group){
  ## This function is designed to obtain the t-value from a logistic regression model, assuming that the 'group' variable is binary (0-1), and 'vdat' is a variable with only one degree of freedom.
  ## The reason for not using t.test to calculate the t-value is that t.test is specifically designed for continuous variables and has distribution requirements.
  ## In contrast, glm can handle both continuous and categorical variables.
  tmpglm = glm(group ~ vdat, family = "binomial")
  glmout = summary(tmpglm)$coefficients[2,3:4]
  return(glmout)
}
