
#' XGpred: Building Risk Classification Predictive Models using Survival Data
#' 
#' @description The XGpred function serves as a wrapper for our novel XGpred algorithm. This algorithm combines XGBoost, clustering, and survival probability prediction
#' to create stable high and low-risk groups based on survival data. These groups are subsequently filtered and leveraged for constructing two models: an XGpred linear prediction score model and an empirical Bayesian-based binary risk classification model.
#' 
#' @details It's important to note that this function does not include a variable selection step. All provided variables will be used in constructing the predictive models.
#' If variable selection is desired, the LASSO_plus function within the csmpv R package can be called.
#'
#' @param data A data matrix or a data frame, samples are in rows, and features/traits are in columns
#' @param varsIn A vector of variables used for prediction model
#' @param time Time variable name 
#' @param event Event variable name
#' @param nrounds An integer to indicate how many times 
#' @param probcut Probability cutoff for risk group classification
#' @param outfile A string for the output file including path if necessary but without file type extension. 
#' @author Aixiang Jiang
#' @return A list is returned with the following seven items:
#' \item{XGBplus}{XGBplus object}
#' \item{stables}{Stabel high and low risk group samples identified by XGBplus}
#' \item{weights}{Weights for each variables used in the model}
#' \item{modelPars}{Mean and standard error of model scores for each risk group}
#' \item{XGpred_score}{Model XGpred score}
#' \item{XGpred_prob}{Empirical Bayesian probability based on model XGpred score}
#' \item{XGpred_prob_class}{Risk group classification based on XGpred_prob for the given probability cutoff}
#' \item{probcut}{Probability cutoff for risk group classification}
#' @references 
#'   Tianqi Chen and Carlos Guestrin, "XGBoost: A Scalable Tree Boosting System", 22nd SIGKDD Conference on Knowledge Discovery and Data Mining, 2016, https://arxiv.org/abs/1603.02754
#'   Aoki T, Jiang A, Xu A et al. The spatially resolved tumor microenvironment predicts treatment outcome in relapsed/refractory Hodgkin lymphoma,
#'          bioRxiv 2023.05.19.541331; doi: https://doi.org/10.1101/2023.05.19.541331
#' @export


XGpred = function(data = NULL, varsIn = NULL, time = NULL, event = NULL, nrounds = 5, probcut = 0.8,
                  outfile = "nameWithPath"){
  
  ######## call XGBplus_cox
  xgbres = XGBplus_cox(data, varsIn, time, event, nrounds)
  
  ## write out the internal validation results for each boosting iteration
  sink(paste0("Internal_validation_",outfile,".txt"))
  xgbres$XGBoost_model$evaluation_log
  sink()
  
  ## stablehighs, stablelows, will get from XGBplus_cox output
  stablehighs = xgbres$stable_class$highs
  stablelows = xgbres$stable_class$lows
  
  stabledat = data.frame(data[c(stablehighs, stablelows),])
  colnames(stabledat) = colnames(data)
  stabledat$class = c(rep(1,length(stablehighs)), rep(0, length(stablelows)))
  #stabledat$class = as.factor(stabledat$class)
  ## the following sapply line only works well for continuous and binary variables in varsIn
  ## when there are more than 2 levels of categorical variables in varsIn, should change them into dummy variables before calling this function
  tout = sapply(stabledat[,varsIn], gettValue, group = stabledat$class) 
  
  pred = data.matrix(data[,colnames(tout)]) %*% tout[1,] 
  
  hpXGpred = pred[stablehighs,1]
  lpXGpred = pred[stablelows,1]
  
  ## change on 20230815,  remove stable classes that have overlapping model scores
  if(min(hpXGpred) < max(lpXGpred)){
    nonsubs = find_non_overlap_sets(lpXGpred, hpXGpred)
    lpXGpred = nonsubs[[1]]
    hpXGpred = nonsubs[[2]]
  }

  gmeans = c(mean(hpXGpred, na.rm = T), mean(lpXGpred, na.rm = T))
  gsds = c(stats::sd(hpXGpred, na.rm = T), stats::sd(lpXGpred, na.rm = T))
  XGpred_prob = getProb(pred, groupMeans = gmeans, groupSds = gsds)  
  
  XGpred_prob_class = ifelse(XGpred_prob >= probcut, "High", "Low")
  ## I need a list, update on 20230815
  stables = list(names(hpXGpred), names(lpXGpred))
  names(stables) = c("stable_high", "stable_low")
  wts = t(tout)
  weights = wts[,1]
  names(weights) = rownames(wts)
  modelPars = cbind(gmeans, gsds)
  colnames(modelPars) = c("mean", "sd")
  rownames(modelPars) = c("stable_high", "stable_low")
  ## pred, XGpred_prob, XGpred_prob_class
  outs = list(xgbres,stables, weights, modelPars, pred, XGpred_prob, XGpred_prob_class, probcut)
  names(outs) = c("XGBplus","stables", "weights", "modelPars", "XGpred_score", "XGpred_prob", "XGpred_prob_class", "probcut")
  return(outs)
}