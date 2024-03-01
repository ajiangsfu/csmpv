#' A function rather aimed at developers
#' @noRd
#' 

find_non_overlap_sets = function(v1, v2){ ## v1 and v2 are two numeric vectors, and max(v1) < max(v2)
  outs = list(v1,v2)
  
  if(max(v1) > min(v2)){
    v1n = length(v1)
    v2n = length(v2)
    vv = sort(c(v1, v2))
    vv1 = vv[1:v1n]
    vv2 = vv[(v1n+1): (v1n+v2n)]
    
    vv1 = subset(vv1, vv1 %in% v1)
    vv2 = subset(vv2, vv2 %in% v2)
    
    outs = list(vv1, vv2)
  }
  
  return(outs)
  
}