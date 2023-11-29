#' A function rather aimed at developers
#' @importFrom stats glm
#' @noRd

gettValue = function(vdat,group){
  ## this function is to get t value fron glm assuming that the group is a 0-1 variable, and vdat is a variable with 1 df only
  ## the reason that I do not use t.test to get t value
  
  tmpglm = glm(group ~ vdat, family = "binomial")
  glmout = summary(tmpglm)$coefficients[2,3:4]
  return(glmout)
}