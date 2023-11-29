#' Predicting Risk Group Classification for a New Data Set
#' 
#' @description The XGpred.predict function is designed to predict risk group classifications for a new data set. This prediction is based on the assumption that an XGpred object is available from a training data set and that the new data set is comparable to the training data set.
#' 
#' In scenarios where the new and training data sets are not directly comparable, a calibration cohort is required to ensure accurate risk group predictions.
#'
#' @param newdat A data matrix or a data frame where samples are in rows and features/traits are in columns. It should include all variables needed
#' for the prediction model.
#' @param XGpredObj An XGpred object returned by the XGpred function. Although not all items are needed for prediction purposes, using the XGpred object as input is convenient.
#' @param scoreShift A calibration value to subtract from the current model score if needed.
#' @author Aixiang Jiang
#' @return A data frame containing XGpred_score, XGpred_prob, and XGpred_prob_class
#' @references 
#'     Aoki T, Jiang A, Xu A et al. The spatially resolved tumor microenvironment predicts treatment outcome in relapsed/refractory Hodgkin lymphoma,
#'          bioRxiv 2023.05.19.541331; doi: https://doi.org/10.1101/2023.05.19.541331
#' @export

XGpred_predict = function(newdat = NULL, XGpredObj = NULL, scoreShift = 0){
  
  wts = XGpredObj$weights
  XGpred_score = data.matrix(newdat[,names(wts)]) %*% wts
  
  XGpred_score = XGpred_score - scoreShift
  
  ## the following line is changed accordingly
  XGpred_prob = getProb(XGpred_score, groupMeans = ( XGpredObj$modelPars)[,1], groupSds = ( XGpredObj$modelPars)[,2])
  
  XGpred_prob_class = ifelse(XGpred_prob >= XGpredObj$probcut, "High", "Low")
  
  ## XGpred_score, XGpred_prob, XGpred_prob_class
  outs = data.frame(cbind(XGpred_score, XGpred_prob))
  colnames(outs) = c("XGpred_score", "XGpred_prob")
  outs$XGpred_prob_class = XGpred_prob_class[,1]
  
  return(outs)
}