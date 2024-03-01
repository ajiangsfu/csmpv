
#' A function rather aimed at developers
#' @importFrom stats dnorm
#' @noRd

getProb = function(inscore, groupMeans, groupSds){
  ### assume groupMeans contain 2 values for 2 group, and the 1st one is for positive group
  ### assume groupSds contain 2 values for 2 group, and the 1st one is for positive group
  #### d1 and d0 are density value for the given observed score under normal assumption
  d1 = dnorm(inscore,mean = groupMeans[1], sd= groupSds[1])  
  d0 = dnorm(inscore,mean = groupMeans[2], sd= groupSds[2])
  return(d1/(d1+d0))

}