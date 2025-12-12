#' XGpred: Building Risk Classification Predictive Models using Survival Data
#' 
#' @description
#' The XGpred function is designed to generate an empirical Bayesian-based binary risk classification model with survival data based on our novel XGpred algorithm, 
#' combining XGBoost and traditional survival analysis.
#' 
#' @details
#' If variable selection is needed, three variable selection methods are provided. Either the given variable or the selected variable list is used to build both an XGBoost model
#' and a traditional Cox model. Risk scores for each model are calculated and ranked, then averaged for each sample. 
#' The top 1/3 of samples are defined as the high-risk group, while the bottom 1/3 of samples are defined as the low-risk group. 
#' The binary risk classification model is built based on these two risk groups using either the given variable or the selected variable list. 
#' The model is a linear combination of these variables, with weights defined as t values derived from the single-variable linear model
#' of each variable on the two groups. Finally, the classification is based on empirical Bayesian probabilities.
#' 

#' @param data A data matrix or a data frame where samples are in rows and features/traits are in columns.
#' @param varsIn A vector of variables used for the prediction model.
#' @param selection Logical. Default is FALSE. If TRUE, three variable selection methods can be chosen.
#' @param vsMethod When "selection" is set to TRUE, three variable selection methods can be chosen, with LASSO2 as the default method. The other two methods
#'                 are "LASSO2plus" and "LASSO_plus."
#' @param time Time variable name.
#' @param event Event variable name.
#' @param nrounds The maximum number of boosting iterations.
#' @param probcut Probability cutoff for risk group classification. Default is set to 0.8.
#' @param nclass Number of risk groups. By default, it is 2; any samples not classified into high-risk groups are classified into the low-risk group. 
#'               When 3 is chosen, samples are classified into low, middle, and high-risk groups.
#' @param topN An integer indicating how many variables to select if LASSO_plus is chosen as the variable selection method.
#' @param outfile A string for the output file, including the path if necessary but without a file type extension.
#' @param nthread The number of parallel threads used to run XGBoost.
#' @param gamma The minimum loss reduction required to make a further partition on a leaf node of the tree.
#' @param max_depth The maximum depth of a tree. Increasing this value will make the model more complex and more likely to overfit.
#' @param outfile A string for the output file, including the path if necessary but without the file type extension.
#' @param eta The step size shrinkage used in the update to prevent overfitting.

#' @import xgboost
#' @import survival
#' @importFrom stats .getXlevels
#' @importFrom stats model.matrix
#' @importFrom stats terms
#' @importFrom stats model.frame
#' 
#' @author Aixiang Jiang
#' @return A list is returned with the following seven items:
#' \item{ranks}{Ranks from XGboost and Cox}
#' \item{twoEnds}{High and low risk group samples identified by mean ranks from XGBoost and Cox models}
#' \item{weights}{Weights for each variables used in the model}
#' \item{modelPars}{Mean and standard error of model scores for each risk group}
#' \item{nclass}{Number of risk groups}
#' \item{XGpred_score}{Model XGpred score}
#' \item{XGpred_prob}{Empirical Bayesian probability based on model XGpred score}
#' \item{XGpred_prob_class}{Risk group classification based on XGpred_prob for the given probability cutoff}
#' \item{probcut}{Probability cutoff for risk group classification}
#' @references 
#'  Tianqi Chen and Carlos Guestrin (2016), "XGBoost: A Scalable Tree Boosting System", 22nd SIGKDD Conference on Knowledge Discovery and Data Mining, 2016, https://arxiv.org/abs/1603.02754
#'  
#'  Aoki T, Jiang A, Xu A et al.,(2023) Spatially Resolved Tumor Microenvironment Predicts Treatment Outcomes in Relapsed/Refractory Hodgkin Lymphoma. J Clin Oncol. 2023 Dec 19:JCO2301115. doi: 10.1200/JCO.23.01115.

#' @examples
#' # Load in data sets:
#' data("datlist", package = "csmpv")
#' tdat = datlist$training
#' 
#' # The function saves files locally. You can define your own temporary directory. 
#' # If not, tempdir() can be used to get the system's temporary directory.
#' temp_dir = tempdir()
#' # As an example, let's define Xvars, which will be used later:
#' Xvars = c("highIPI", "B.Symptoms", "MYC.IHC", "BCL2.IHC", "CD10.IHC", "BCL6.IHC")

#' # For given time-to-event outcome and Xvars, we can build up a binary risk classification:
#'  xgobj = XGpred(data = tdat, varsIn = Xvars, 
#'                 time = "FFP..Years.", event = "Code.FFP", 
#'                 outfile = paste0(temp_dir, "/XGpred"))
#' # You might save the files to the directory you want.
#' 
#' # To delete the temp_dir, use the following:
#' unlink(temp_dir)

#' @export

XGpred = function(data = NULL, varsIn = NULL, selection = FALSE, 
                  vsMethod = c("LASSO2", "LASSO2plus", "LASSO_plus"),
                  time = NULL, event = NULL, nrounds = 5, probcut = 0.8,
                  nthread = 2, gamma = 1, max_depth = 3, eta = 0.3,
                  nclass = c(2,3), topN = 10, outfile = "nameWithPath"){
  nclass = nclass[1]
  ## variable selection
  if(selection){
    vsMethod = vsMethod[1]
    if(vsMethod == "LASSO_plus"){
      afit = LASSO_plus(data = data,  biomks = varsIn, outcomeType = "time-to-event", 
                        time = time, event = event, outfile = outfile, topN = topN)
      aform = afit$fit$formula
      tmp = strsplit(toString(aform), split = ",")[[1]]
      tmp = tmp[length(tmp)]
      tmp = strsplit(tmp, split = "\\+")[[1]]
      if(length(tmp) > 1){
        tmp = sapply(tmp, function(xx){
          gsub(" ", "", xx)
        })
      }else{
        tmp = gsub(" ", "", tmp)
      }
      varsIn = tmp
    }else if(vsMethod == "LASSO2plus"){
      afit = LASSO2plus(data = data,  biomks = varsIn, outcomeType = "time-to-event", 
                        time = time, event = event, outfile = outfile)
      aform = afit$fit$formula
      tmp = strsplit(toString(aform), split = ",")[[1]]
      tmp = tmp[length(tmp)]
      tmp = strsplit(tmp, split = "\\+")[[1]]
      if(length(tmp) > 1){
        tmp = sapply(tmp, function(xx){
          gsub(" ", "", xx)
        })
      }else{
        tmp = gsub(" ", "", tmp)
      }
      varsIn = tmp
    }else{ ## LASSO2 as the default method
      lasres = LASSO2(data = data, biomks = varsIn, outcomeType = "time-to-event", 
                          time = time, event = event, outfile = outfile)
      varsIn = names(lasres$coefs)
    }
  }
  
  x.train <- data[, c(varsIn, time, event), drop = FALSE]
  
  # Create formula from biomarkers
  formula <- as.formula(paste("~", paste(varsIn, collapse = " + ")))
  
  # Create model matrix - this properly handles factors
  tdat <- model.matrix(formula, data = x.train)
  
  # Remove intercept (if present) and convert to matrix
  tdat <- tdat[, -1, drop = FALSE]  # Remove intercept column
  storage.mode(tdat) <- "double"
  
  # Store factor information for consistent prediction later
  factor_info <- list(
    xlevels = .getXlevels(terms(formula), x.train),
    contrasts = attr(tdat, "contrasts")
  )
  
  # Survival label: time for events, -time for censored
  time_vals <- x.train[[time]]
  event_vals <- x.train[[event]]
  survival_label <- time_vals * (2 * as.numeric(event_vals) - 1)
  
  # Create DMatrix - CRAN-safe approach
  # Option 1: Direct approach (usually works)
  Dtrain <- xgboost::xgb.DMatrix(data = tdat, label = survival_label)


  # Parameters for xgb.train()
  params <- list(
    objective = "survival:cox",
    nthread = nthread,
    gamma = gamma,
    max_depth = max_depth,
    eta = eta
  )
  
  modeln = xgboost::xgb.train(
    params = params,
    data = Dtrain,
    nrounds = nrounds,
    verbose = 2
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
  ranks = data.frame(ranks)
  
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
  
  if(nclass == 3){
    XGpred_prob_class = ifelse(XGpred_prob >= probcut, "High", ifelse(XGpred_prob <= 1-probcut, "Low", "Middle"))
  }else{
    XGpred_prob_class = ifelse(XGpred_prob >= probcut, "High", "Low")
  }
  
  ## I need a list
  twoEnds = list(names(hpXGpred), names(lpXGpred))
  names(twoEnds) = c("highs", "lows")
  wts = t(tout)
  weights = wts[,1]
  names(weights) = rownames(wts)
  modelPars = cbind(gmeans, gsds)
  colnames(modelPars) = c("mean", "sd")
  rownames(modelPars) = c("highs", "lows")
  ## pred, XGpred_prob, XGpred_prob_class
  outs = list(ranks,twoEnds, weights, modelPars,nclass, pred, XGpred_prob, XGpred_prob_class, probcut)
  names(outs) = c("ranks","twoEnds", "weights", "modelPars","nclass", "XGpred_score", "XGpred_prob", "XGpred_prob_class", "probcut")
  return(outs)
}