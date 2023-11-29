
#' A function rather aimed at developers
#' @noRd

standardize = function(dataIn, byrow = TRUE){
  ind = ifelse(byrow == TRUE, 1,2)
  outs = apply(dataIn, ind, FUN = function(xx){
    mm = mean(xx, na.rm = TRUE)
    ss = stats::sd(xx, na.rm = TRUE)
    (xx-mm)/ss 
  })
  
  if(byrow == TRUE){
    outs = t(outs)
  }
  return(outs)

}