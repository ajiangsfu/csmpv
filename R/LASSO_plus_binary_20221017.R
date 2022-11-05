
#' LASSO_plus on variable selection for binary outcome
#' @description This is the function to select variables that are associated with a bianry outcome variable based on my new algorithm: LASSO_plus. 
#' @details This function is to use LASSO_plus algorithm to select variables that are assoicated with on an binary outcome variable for a given data set.
#' LASSO_plus selects variables based on LASSO, single variable regession, and stepwise regression.  
#' @param data A data matrix or a data frame, samples are in columns, and features/traits are in rows.
#' @param biomks A vector of potential biomarkers for variable selection, they should be a subset of "data" column names. 
#' @param Y Outcome varialbe name when the outcome type is either "binary" or "continuous". 
#' @param time Time variable name when outcome type is "time-to-event".
#' @param event Event variable name when outcome type is "time-to-event".
#' @keywords variable selection
#' @author Aixiang Jiang
#' @return A list with three items is returned:
#' 
#' do it later
#' 
#' \item{name}{content}
#' \item{}{}
#' @references 
#' 
#' @export
#' 
#' wrote code first, later, update the above section
LASSO_plus_binary = function(data, Y, biomks, filename = "default", 
                                     title = "default", topN = 10,  forestHeight = 4){
  pid = dat[, patientID]
  morevar = c(y01)
  fdat = dat[,morevar]
  dat = dat[,biomks]
  
  # ### remove variables if their sd < 0.0000001
  sdcol = apply(dat,2,sd, na.rm = TRUE)
  sdtmp = which(sdcol>0.0000001)
  dat = dat[,sdtmp]
  
  vars = intersect(biomks, colnames(dat))
  dat = cbind(fdat, dat)
  
  ### remove NA for Y side if there are any
  natmp = which(is.na(dat$fdat))
  
  if(length(natmp) > 0){
    dat = dat[-natmp,]
  }

  lassoF = glmnet(x= data.matrix(dat[,vars]), y=as.integer(dat$fdat), family = "binomial")
  
  allcoefs = data.matrix(lassoF$beta)
  
  ######################################################
  ### more output for LASSO on 0827
  outfile = paste0("glm_",topN,"_N_",filename, "_", title)
  lassoPlot = paste0(outfile, "_LASSO.pdf")
  lassoCoef = gsub("pdf", "csv",lassoPlot)
  
  pdf(lassoPlot)
  plot(lassoF)
  #plotres(lassoF)
  plot(lassoF,xvar="lambda",label=TRUE)
  dev.off()
  
  ####################################################

  ## calculate how many variables are remained, and a0 is in a different output 
  coef01 = apply(allcoefs, c(1,2), FUN = function(xx) { ifelse(abs(xx) > 0.00001, 1, 0)})
  csum = colSums(coef01)
  
  tt = csum[duplicated(csum)]
  
  ### to find who close to top N
  tmp = which.min(abs(tt - topN))  ### it gives the 1st matching
  
  ### now, find lambda for this tmp
  scoefs = allcoefs[,names(tmp)]
  
  #############################
  ## write out all coefs and lambda as well on 0827
  colnames(allcoefs) = paste(colnames(allcoefs), lassoF$lambda, sep="_lambda")
  write.csv(allcoefs, lassoCoef)
  
  ### write out scoefs, 20220204
  write.csv(scoefs, paste0("LASSOselected_",lassoCoef))
  
  ##############################
  ### get out the variable names whose coefs are not 0
  topnsh = which(abs(scoefs) > 0.00001)
  
  ### get variable names only
  topnsh = names(topnsh)  ### this is the variable list from LASSO
  
  if(length(topnsh) < 1){
    stop("No significant varaible was selected by LASSO")
  }
  
  f1 =  as.formula(paste("fdat", paste(topnsh, collapse=" + "), sep=" ~ "))
  
  fit1 = glm(formula = f1, family="binomial", data = dat)  
  
  out1 = summary(fit1)$coef
  out2 = exp(cbind("Odds_ratio" = coef(fit1), confint.default(fit1, level = 0.95)))
  
  tmp = intersect(rownames(out1), rownames(out2))
  
  if(length(tmp) > 1){
    fitout = cbind(out1[tmp,], out2[tmp,])
  }else{
    stop("No significant varaible was selected by LASSO")
  }
  
  ### remove intercept, order by odds ratio
  fitout = fitout[-1,]
  
  fitout = fitout[order(fitout[,5], decreasing = TRUE),]
  
  outfile = paste0(topN,"_N_",y01,"_",filename, "_", title)
  outfile1 = paste0(outfile,"logodds.pdf")
  
  pdf(outfile1, width = 13, height = forestHeight)
  log10forestPlot(resultIn = fitout, effInds = 5:7, breaks = c(0.1, 0.5, 1, 2, 10))
  dev.off()
  
  outfile2 = paste0(outfile,"_glm.csv")
  write.csv(fitout, outfile2)  
  
  ### model score
  modelScore = predict(fit1, newdata = dat, type="response")
  ### should add modelScore with other info
  modelScore = cbind(fdat, pid, modelScore)
  colnames(modelScore) = c(y01, patientID, "Model_score")
  write.csv(modelScore, gsub("glm", "modelScore", outfile2))
  
  
  ###############################################################################
  ###################### the above are customized LASSO results ################
  ###############################################################################
  
  ### now, let's try to use single variable approach
  #survY = paste0("Surv(", morevar[1],",", morevar[2], ")")  ## this is duplicated, and should remove
  bout = t(sapply(vars, function(var1){
    out1 = glm(as.formula(paste("fdat", var1, sep=" ~ ")), family="binomial", data = dat)
    return(summary(out1)$coefficients[2,])
  }))
  
  colnames(bout) = c("coef", "se(coef)", "z", "p")
  ### get FDR
  bout = data.frame(bout)
  bout$adjusteP = p.adjust(bout$p, method = "BH")
  
  ### order them and write out
  bout = bout[order(bout$adjusteP),]
  ## cutoff, FDR < 0.05
  bout = subset(bout, bout$adjusteP <= 0.05)
  if(dim(bout)[1] > 0){ 
    bout = head(bout, n = topN)
    write.csv(bout, gsub("LASSO", "singleVar", lassoCoef))
  }  
  topn2 = rownames(bout)  ### even this is empty, the following union still works
  
  ### combine topn1 and topn2, run a full model and step
  topn1 = union(topnsh,topn2)
  
  if(length(topn1) < 1){
    stop("No significant varaible was selected by LASSO_plus")
  }
  
  vX = paste(topn1, collapse=" + ")
  
  fit1 = glm(as.formula(paste("fdat", vX, sep=" ~ ")), family="binomial",data = dat)
  
  ###########################################################
  ####### the following are LASSO_plus outputs ################
  
  fit2 = step(fit1)
  
  afit1 = summary(fit2)$coef
  ## I need odds ratio as well
  afit2 = exp(cbind("Odds_ratio" = coef(fit2), confint.default(fit2, level = 0.95)))
  
  tmp = intersect(rownames(afit1), rownames(afit2))
  
  if(length(tmp) > 1){
    fitout = cbind(afit1[tmp,], afit2[tmp,])
  }else{
    stop("No significant varaible was selected by LASSO_plus")
  }
  
  ### remove intercept, order by odds ratio
  fitout = fitout[-1,]
  
  fitout = fitout[order(fitout[,5], decreasing = TRUE),]
  
  outfile1 = paste0(outfile,"LASSO_plus.csv")
  outfile2 = paste0(outfile,"LASSO_plus_forest.pdf")
  
  write.csv(fitout, outfile1)
  
  pdf(outfile2, width = 13, height = forestHeight)
  log10forestPlot(resultIn = fitout, HRvsORplot = FALSE, effInds = 5:7, breaks = c(0.1, 0.5, 1, 2, 10))
  dev.off()

  ### re-run fit1 based on variables in fit2
  vX = paste(rownames(fitout), collapse=" + ")
  fit1 =  glm(as.formula(paste("fdat", vX, sep=" ~ ")), family="binomial",data = dat)
  
  ## model checking
  fpdf = gsub("forest", "ModelCheck", outfile2)
  pdf(fpdf, width = 8, height = 8)
  plotGlm(fit1)
  dev.off()
  
  ### return
  return(vif(fit1))
  
}


### I put the R code file and the data file under the same folder:
library(rstudioapi)
current_working_dir = dirname(rstudioapi::getActiveDocumentContext()$path)

setwd(current_working_dir)

load("glmdattest.RData")
# > ls()
# [1] "binaryVarSelectionLASSO_plus" "codetest"                    
# [3] "coxphOutput"                  "current_working_dir"         
# [5] "glmdat"                       "log10forestPlot"             
# [7] "plotGlm"                      "plotTable"                   
# [9] "variabletest" 

binaryVarSelectionLASSO_plus(as.data.frame(glmdat), codetest, variabletest, "res_id", "test", "glmLASSOoutput")



# 
# ### test code
# dat = as.data.frame(glmdat)
# y01 = codetest
# biomks = variabletest
# filename = "test"
# title = "glmLASSOoutput"


####### old code  ##########
# tdat = read.csv("topMean_ASCT1_postBMTFFS_49ids_20220503.csv",
#                 header = T, row.names = 1, stringsAsFactors = F)
# vars = grep("_spatial_score", colnames(tdat))
# vars = colnames(tdat)[vars]
# 
# ## in the following: y01 = "CODE_postBMTFFS" is not a good example, but I just want to show how the function works
# res = binaryVarSelectionLASSO_plus(dat=tdat, y01 = "CODE_postBMTFFS", biomks = vars,topN = 10,
#                                 filename = "topMean", title = "49ids", forestHeight = 3)
# 
# # > res  ## VIF, there are some VIF > 10, indicates potential multicollinearity variables
# # CXCR5_HRS_spatial_score PDL1CD80_Mac_Myel_spatial_score          PDL1_CD4_spatial_score               Mac_spatial_score 
# # 2.102444                        4.764440                        3.428928                        9.723370 
# # LAG3_CD4_spatial_score              Endo_spatial_score           CXCR5_B_spatial_score         TIM3_Treg_spatial_score 
# # 10.402737                        3.971172                        3.414042                       14.561848 
# # TIM3_Mac_spatial_score 
# # 2.351734 
# 
# 
# 
# 
# 
# 
# 
# 
