#' A function rather aimed at developers
#' @import mclust 
#' @import survival
#' @import xgboost
#' 
#' @noRd

XGBplus_cox = function(data, varsIn, time, event, nrounds = 5 ) {
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
  
  pn = stats::predict(modeln, Dtrain)
  
  ### next step, find stable highs and lows with 3 methods, should put all of them into the above R function, with given XGBoost score, 
  ## 1) mclust
  #clures1 = mclust::Mclust(pn, G = 3)
  clures2 = mclust::Mclust(log(pn), G = 3)
  
  res = data.frame(cbind(data[,c(time, event)], pn))
  colnames(res)[3] = "XGBoostScore"
  
  #res$XGBoost_Mclust = clures1$classification  ## changed on 20230814, I do not need this any more
  #res$logXGBoost_Mclust = clures2$classification
  #res$XGBoost_Mclust_class = ifelse(res$logXGBoost_Mclust == res$XGBoost_Mclust, res$XGBoost_Mclust, NA) 
  # big change on 20230814: 
  # only keep result from logXGBoost_Mclust: this is because xgboost predict for time to event is exp(linear model score), so we should cluster the model score after log transformation
  res$XGBoost_Mclust_class = clures2$classification
  res$XGBoost_Mclust_class = factor(res$XGBoost_Mclust_class)
  
  ### XGBoost_Mclust13 should based on two ends of XGBoost_Mclust_class
  c1 = subset(res, res$XGBoost_Mclust_class == 1)
  c2 = subset(res, res$XGBoost_Mclust_class == 2)
  c3 = subset(res, res$XGBoost_Mclust_class == 3)
  rc1 = range(c1$XGBoostScore)
  rc2 = range(c2$XGBoostScore)
  rc3 = range(c3$XGBoostScore)
  
  rc = rbind(rc1,rc2,rc3)
  max1 = min(rc[,2])
  min3 = max(rc[,1])
  res$XGBoost_Mclust13 = ifelse(res$XGBoostScore <= max1, 0, ifelse(res$XGBoostScore >= min3, 1, NA))
  res$XGBoost_Mclust13 = factor(res$XGBoost_Mclust13)
  
  ## 2) change on 20230629
  ## calculate survival prob based on exp(linear model score), directly from XGBoost
  coxfit = coxph(Surv(data[,time], data[,event]) ~ 1)
  bh = basehaz(coxfit, centered=FALSE)
  ## I need time for baseline hazard, use two years for now?
  tmp = which.min(abs(bh[,2]-2)) ## this is not exactly, but roughly correct
  bh = bh[tmp,1]
  
  ## # Calculate the survival probability using the baseline hazard and XGBoost prediction
  survp = exp(-bh*pn) ## pn = exp(linear prediction score), from XGBoost
  
  res$xgb.survProb = survp
  
  lowcut =  max1 ## based on above Mclust results
  
  ## change on 20230815
  highcut = min(2, min3) ## 2 from experience based on surv curve based on RHL4S
  res$XGBoost_Mclust_surv = ifelse(res$XGBoostScore >= highcut, 1, ifelse(res$XGBoostScore <= lowcut, 0, NA))
  
  ## extra step:
  ## also, 0.5 is a natural cutoff, to be conservative, use 0.3 as high end cutoff for samples that are not in high risk group
  res$XGBoost_Mclust_surv = ifelse(is.na(res$XGBoost_Mclust_surv) & res$xgb.survProb < 0.3, 1, res$XGBoost_Mclust_surv) 

  highs = rownames(subset(res, res$XGBoost_Mclust_surv == 1))
  lows = rownames(subset(res, res$XGBoost_Mclust_surv == 0))
  
  ## 3) survival data: time, notice that here time is a variable to indicate actually column name for time
  tmp = subset(res[lows,], res[lows,time] >=5)
  
  ## define stable highs and lows
  lows = rownames(tmp)
  
  stables = list(highs, lows)
  names(stables) = c("highs", "lows")
  outs = list(modeln, res, stables)
  names(outs) = c("XGBoost_model", "XGplus_cox_result", "stable_class")
  return(outs)
}