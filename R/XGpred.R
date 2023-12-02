#' XGpred: Building Risk Classification Predictive Models using Survival Data
#' 
#' @description
#' The XGpred function is designed to generate an empirical Bayesian-based binary risk classification model with survival data based on our novel XGpred algorithm, 
#' combining XGBoost and traditional survival analysis.
#' 
#' @details
#' If variable selection is needed, the LASSO2 algorithm is processed. The given variable or the LASSO2 selected variable list is used to build both an XGBoost model
#' and a traditional Cox model. Risk scores for each model are calculated and ranked separately. The ranks are then averaged for each sample. 
#' The top 1/3 samples are defined as the high-risk group, while the low 1/3 samples are defined as the low-risk group. 
#' The binary risk classification model is built based on the two risk groups using the given variable or the LASSO2 selected variable list. 
#' The model is a linear combination of these variables, with weights as t values from the single-variable linear model of each variable on the two groups. 
#' Finally, the classification is based on empirical Bayesian probabilities.
#'#' 
#' @param data A data matrix or a data frame, samples are in rows, and features/traits are in columns
#' @param varsIn A vector of variables used for prediction model
#' @param selection Logical. If TRUE, variable selection is performed using LASSO2.
#' @param time Time variable name 
#' @param event Event variable name
#' @param nrounds An integer to indicate how many times 
#' @param probcut Probability cutoff for risk group classification
#' @param outfile A string for the output file including path if necessary but without file type extension. 
#' @import xgboost
#' @import survival
#' @author Aixiang Jiang
#' @return A list is returned with the following seven items:
#' \item{ranks}{Ranks from XGboost and Cox}
#' \item{twoEnds}{High and low risk group samples identified by mean ranks from XGBoost and Cox models}
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


XGpred = function(data = NULL, varsIn = NULL, selection = FALSE, time = NULL, event = NULL, nrounds = 5, probcut = 0.8,
                  outfile = "nameWithPath"){
  
  ## variable selection
  if(selection){
    
    lasres = LASSO2(data = data, biomks = varsIn, outcomeType = "time-to-event", 
                    time = time, event = event, outfile = outfile)
    varsIn = names(lasres$coefs)
  }
  
  x.train = data[,c(varsIn, time, event)]
  ### XGB
  num_feature = dim(x.train)[2]
  x.train.xgb = data.matrix(x.train)
  dtrain = list(data=x.train.xgb[,c(1:(num_feature-2))],label=x.train.xgb[,(num_feature-1)]*(-(-1)^(as.numeric(x.train.xgb[,num_feature]))))	
  Dtrain = xgboost::xgb.DMatrix(dtrain$data,label=dtrain$label)
  modeln = xgboost::xgboost(
    objective = "survival:cox",
    data = Dtrain,
    nrounds = nrounds,
    nthread = 2, ## this is not important for small data set
    verbose = 2 
    # If 0, xgboost will stay silent. If 1, xgboost will print information of performance. 
    # If 2, xgboost will print information of both performance and construction progress information
  )
  
  pn = stats::predict(modeln, Dtrain)  ## this is risk score: exp(linear prediction)
  xranks = rank(pn)
  
  ## run a cox model as well
  survY = paste0("survival::Surv(", time,",", event, ")")
  survX = paste(varsIn, collapse=" + ")
  coxres = survival::coxph(as.formula(paste(survY, survX, sep=" ~ ")), data = data)
  
  ## the risk score exp(lp) ("risk"), same format as in pn (xgboost predicted values)
  cscores = stats::predict(coxres, type="risk")
  cranks = rank(cscores)
  
  ranks = cbind(xranks, cranks)
  rownames(ranks) = rownames(data)
  rmeans = rowMeans(ranks)
  ranks = cbind(ranks, rmeans)
  
  ## get quantiles: 1/3 and 2/3
  q1 = stats::quantile(ranks$rmeans,1/3)
  q3 = stats::quantile(ranks$rmeans,2/3)
  
  highs = rownames(subset(ranks, ranks$rmeans >= q3))
  lows = rownames(subset(ranks, ranks$rmeans <= q1))
  
  colnames(ranks) = c("rank_XGBoost", "ranks_Cox", "mean")
  
  sdat = data.frame(data[c(highs, lows),])
  colnames(sdat) = colnames(data)
  sdat$class = c(rep(1,length(highs)), rep(0, length(lows)))
  tout = sapply(sdat[,varsIn], gettValue, group = sdat$class) 
  
  pred = data.matrix(data[,colnames(tout)]) %*% tout[1,] 
  
  hpXGpred = pred[highs,1]
  lpXGpred = pred[lows,1]
  
  gmeans = c(mean(hpXGpred, na.rm = T), mean(lpXGpred, na.rm = T))
  gsds = c(stats::sd(hpXGpred, na.rm = T), stats::sd(lpXGpred, na.rm = T))
  XGpred_prob = getProb(pred, groupMeans = gmeans, groupSds = gsds)  
  
  XGpred_prob_class = ifelse(XGpred_prob >= probcut, "High", "Low")
  ## I need a list, update on 20230815
  twoEnds = list(names(hpXGpred), names(lpXGpred))
  names(twoEnds) = c("highs", "lows")
  wts = t(tout)
  weights = wts[,1]
  names(weights) = rownames(wts)
  modelPars = cbind(gmeans, gsds)
  colnames(modelPars) = c("mean", "sd")
  rownames(modelPars) = c("highs", "lows")
  ## pred, XGpred_prob, XGpred_prob_class
  outs = list(ranks,twoEnds, weights, modelPars, pred, XGpred_prob, XGpred_prob_class, probcut)
  names(outs) = c("ranks","twoEnds", "weights", "modelPars", "XGpred_score", "XGpred_prob", "XGpred_prob_class", "probcut")
  return(outs)
}