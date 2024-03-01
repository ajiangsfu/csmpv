#' All-in-one Modelling with csmpv R package
#'
#' @description
#' This function is designed to simplify the process of building, evaluating and comparing different modelling methods. It offers the flexibility to perform one 
#' or all of the following modelling methods: LASSO2, LASSO2 + regression, LASSO_plus, LASSO2plus, XGBoost, LASSO2 + XGBoost, LASSO_plus + XGBoost, 
#' and LASSO2plus + XGBoost. The models are trained on the training data, and their performance is validated on a separate validation dataset.
#'
#' @details
#' By default, this function runs all eight different modeling methods. However, users can specify the "methods" parameter to choose and run a specific modelling 
#' method of their choice. For clarity, when providing a 'vdat' argument, the function assumes that it contains the outcome variable, and it proceeds with model validation.
#' 
#' @param tdat Training data. It can not be null.
#' @param vdat Validation data. It should contain the same variables as in the training data, including outcome variables. No validation result is saved if it is NULL.
#' @param Ybinary Binary outcome variable for classification.
#' @param varsBinary Names of binary predictors.
#' @param Ycont Continuous outcome variable for regression.
#' @param varsCont Names of continuous predictors.
#' @param time Time-to-event variable for survival analysis.
#' @param event Event/censoring indicator for survival analysis.
#' @param varsSurvival Names of predictors for survival analysis.
#' @param methods Method(s) to use for modeling. If "all," models for all eight methods will be built. 
#'                Otherwise, provide one of the following method names:
#'                - "LASSO2": Variable selection using LASSO2 with a minimum of two remaining variables.
#'                - "LASSO2_reg": Variables selected from LASSO2, followed by regular regression.
#'                - "LASSO_plus": Variables selected from LASSO_plus, followed by regular regression.
#'                - "LASSO2plus": Variables selected from LASSO2plus, followed by regular regression.
#'                - "XGBoost": XGBoost model built without variable selection.
#'                - "LASSO2_XGBoost": Variables selected from LASSO2, followed by XGBoost.
#'                - "LASSO_plus_XGBoost": Variables selected from LASSO_plus, followed by XGBoost.
#'                - "LASSO2plus_XGBoost": Variables selected from LASSO2plus, followed by XGBoost.
#' @param outfileName Prefix for output file names.
#' @author Aixiang Jiang
#' @return A list of trained models and prediction objects. Results are saved to local files.
#' 
#' @examples
#' # Load in data sets:
#' data("datlist", package = "csmpv")
#' tdat = datlist$training
#' vdat = datlist$validation
#' 
#' # The confirmVars function saves files locally. You can define your own temporary directory. 
#' # If not, tempdir() can be used to get the system's temporary directory.
#' temp_dir = tempdir()
#' 
#' # As an example, let's define Xvars, which will be used later:
#' Xvars = c("highIPI", "B.Symptoms", "MYC.IHC", "BCL2.IHC", "CD10.IHC", "BCL6.IHC")

#' # The default setting of this single function generates all models and provides predictions
#' # and validations for each of them. 
#' # Of course, we can also use this all-in-one function to work on one outcome type 
#' # and one model at a time, for example:
#' DZlassoreg = csmpvModelling(tdat = tdat, vdat = vdat,
#'                            Ybinary = "DZsig", varsBinary = Xvars,
#'                            methods = "LASSO2_reg",
#'                            outfileName= paste0(temp_dir, "/just_one"))
#' # This is equivalent to using LASSO2_reg for modeling, followed by prediction and validation 
#' # with rms_model for the classification task "DZsig".
#' # Six result files are then saved locally.
#' # You might want to save the files to the directory you prefer.
#' 
#' # To delete the "temp_dir", use the following:
#' unlink(temp_dir)

#' @export

csmpvModelling = function(tdat=NULL, vdat=NULL, Ybinary =NULL, varsBinary=NULL, Ycont=NULL, varsCont=NULL, time=NULL, event=NULL, varsSurvival=NULL, 
                          methods = c("all", "LASSO2", "LASSO2_reg", "LASSO_plus", "LASSO2plus",
                                      "XGBoost", "LASSO2_XGBoost", "LASSO_plus_XGBoost", "LASSO2plus_XGBoost"), 
                          outfileName = NULL){
 
  outnames = c("bl","cl","tl","blr","clr","tlr","bfit", "cfit","tfit","b2fit", "c2fit", "t2fit",
               "bxfit", "cxfit", "txfit", "blxfit", "clxfit", "tlxfit", "blpxfit", "clpxfit", "tlpxfit",  "bl2xfit", "cl2xfit", "tl2xfit")
  outnames = c(outnames, paste0("p", outnames))
  
  for (name in outnames) {
    assign(name, NULL)
  }
  
  methods = methods[1]
  
  ## binary
  if(!is.null(Ybinary)){
    if (methods == "LASSO2"){
      bl = LASSO2(data = tdat, biomks = varsBinary, Y = Ybinary, outfile = paste0(outfileName,"_binary_LASSO2"))
      if(!is.null(vdat)){
        pbl = LASSO2_predict(bl, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_binary_LASSO2_validate"))
      }
    }else if (methods == "LASSO2_reg"){
      blr = LASSO2_reg(data = tdat, biomks = varsBinary, Y = Ybinary, outfile = paste0(outfileName,"_binary_LASSO2reg"))
      if(!is.null(vdat)){
        pblr = rms_model(blr$fit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_binary_LASSO2reg_validate"))
      }
    }else if (methods == "LASSO_plus"){
      bfit = LASSO_plus(data = tdat, biomks = varsBinary, Y = Ybinary, outfile = paste0(outfileName,"_binary_LASSOplus"))
      if(!is.null(vdat)){
        pbfit = rms_model(bfit$fit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_binary_LASSOplus_validate"))
      }
    }else if (methods == "LASSO2plus"){
      b2fit = LASSO2plus(data = tdat, biomks = varsBinary, Y = Ybinary, outfile = paste0(outfileName,"_binary_LASSO2plus"))
      if(!is.null(vdat)){
        pb2fit = rms_model(b2fit$fit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_binary_LASSO2plus_validate"))
      }
    }else if (methods == "XGBoost"){
      bxfit = XGBtraining(data = tdat, biomks = varsBinary, Y = Ybinary, outfile = paste0(outfileName,"_binary_XGBoost"))
      if(!is.null(vdat)){
        pbxfit = XGBtraining_predict(bxfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_binary_XGBoost_validate"))
      }
    }else if (methods == "LASSO2_XGBoost"){
      blxfit = LASSO2_XGBtraining(data = tdat, biomks = varsBinary, Y = Ybinary, outfile = paste0(outfileName,"_binary_LASSO2_XGBoost"))
      if(!is.null(vdat)){
        pblxfit = XGBtraining_predict(blxfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_binary_LASSO2_XGBoost_validate"))
      }
    }else if (methods == "LASSO_plus_XGBoost"){
      blpxfit = LASSO_plus_XGBtraining(data = tdat, biomks = varsBinary, Y = Ybinary, outfile = paste0(outfileName,"_binary_LASSOplus_XGBoost"))
      if(!is.null(vdat)){
        pblpxfit = XGBtraining_predict(blpxfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_binary_LASSOplus_XGBoost_validate"))
      }
    }else if(methods == "LASSO2plus_XGBoost"){
      bl2xfit = LASSO2plus_XGBtraining(data = tdat, biomks = varsBinary, Y = Ybinary, outfile = paste0(outfileName,"_binary_LASSO2plus_XGBoost"))
      if(!is.null(vdat)){
        pbl2xfit = rms_model(bl2xfit$fit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_binary_LASSO2plus_validate"))
      }
    }else{
      bl = LASSO2(data = tdat, biomks = varsBinary, Y = Ybinary, outfile = paste0(outfileName,"_binary_LASSO2"))
      blr = LASSO2_reg(data = tdat, biomks = varsBinary, Y = Ybinary, outfile = paste0(outfileName,"_binary_LASSO2reg"))
      bfit = LASSO_plus(data = tdat, biomks = varsBinary, Y = Ybinary, outfile = paste0(outfileName,"_binary_LASSOplus"))
      b2fit = LASSO2plus(data = tdat, biomks = varsBinary, Y = Ybinary, outfile = paste0(outfileName,"_binary_LASSO2plus"))
      bxfit = XGBtraining(data = tdat, biomks = varsBinary, Y = Ybinary, outfile = paste0(outfileName,"_binary_XGBoost"))
      blxfit = LASSO2_XGBtraining(data = tdat, biomks = varsBinary, Y = Ybinary, outfile = paste0(outfileName,"_binary_LASSO2_XGBoost"))
      blpxfit = LASSO_plus_XGBtraining(data = tdat, biomks = varsBinary, Y = Ybinary, outfile = paste0(outfileName,"_binary_LASSOplus_XGBoost"))
      bl2xfit = LASSO2plus_XGBtraining(data = tdat, biomks = varsBinary, Y = Ybinary, outfile = paste0(outfileName,"_binary_LASSO2plus_XGBoost"))
      if(!is.null(vdat)){
        pbl = LASSO2_predict(bl, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_binary_LASSO2_validate"))
        pblr = rms_model(blr$fit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_binary_LASSO2reg_validate"))
        pbfit = rms_model(bfit$fit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_binary_LASSOplus_validate"))
        pb2fit = rms_model(b2fit$fit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_binary_LASSO2plus_validate"))
        pbxfit = XGBtraining_predict(bxfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_binary_XGBoost_validate"))
        pblxfit = XGBtraining_predict(blxfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_binary_LASSO2_XGBoost_validate"))
        pblpxfit = XGBtraining_predict(blpxfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_binary_LASSOplus_XGBoost_validate"))
        pbl2xfit = XGBtraining_predict(bl2xfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_binary_LASSO2plus_XGBoost_validate"))
      }
    }
  }
  ## continuous
  if(!is.null(Ycont)){
    if (methods == "LASSO2"){
      cl = LASSO2(data = tdat, biomks = varsCont, outcomeType = "continuous", Y = Ycont, outfile = paste0(outfileName,"_continuous_LASSO2"))
      if(!is.null(vdat)){
        pcl = LASSO2_predict(cl, newdata = vdat, newY = TRUE, outfile =  paste0(outfileName,"_continuous_LASSO2_validate"))
      }
    }else if (methods == "LASSO2_reg"){
      clr = LASSO2_reg(data = tdat, biomks = varsCont, outcomeType = "continuous", Y = Ycont, outfile = paste0(outfileName,"_continuous_LASSO2reg"))
      if(!is.null(vdat)){
        pclr = rms_model(clr$fit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_continuous_LASSO2reg_validate"))
      }
    }else if (methods == "LASSO_plus"){
      cfit = LASSO_plus(data = tdat, biomks = varsCont, outcomeType = "continuous", Y = Ycont, outfile = paste0(outfileName,"_continuous_LASSOplus"))
      if(!is.null(vdat)){
        pcfit = rms_model(cfit$fit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_continuous_LASSOplus_validate"))
      }
    }else if (methods == "LASSO2plus"){
      c2fit = LASSO2plus(data = tdat, biomks = varsCont, outcomeType = "continuous", Y = Ycont, outfile = paste0(outfileName,"_continuous_LASSO2plus"))
      if(!is.null(vdat)){
        pc2fit = rms_model(c2fit$fit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_continuous_LASSO2plus_validate"))
      }
    }else if (methods == "XGBoost"){
      cxfit = XGBtraining(data = tdat, biomks = varsCont, outcomeType = "continuous", Y = Ycont, outfile = paste0(outfileName,"_continuous_XGBoost"))
      if(!is.null(vdat)){
        pcxfit = XGBtraining_predict(cxfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_continuous_XGBoost_validate"))
      }
    }else if (methods == "LASSO2_XGBoost"){
      clxfit = LASSO2_XGBtraining(data = tdat, biomks = varsCont, outcomeType = "continuous", Y = Ycont, outfile = paste0(outfileName,"_continuous_LASSO2_XGBoost"))
      if(!is.null(vdat)){
        pclxfit = XGBtraining_predict(clxfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_continuous_LASSO2_XGBoost_validate"))
      }
    }else if (methods == "LASSO_plus_XGBoost"){
      clpxfit = LASSO_plus_XGBtraining(data = tdat, biomks = varsCont, outcomeType = "continuous", Y = Ycont, outfile = paste0(outfileName,"_continuous_LASSOplus_XGBoost"))
      if(!is.null(vdat)){
        pclpxfit = XGBtraining_predict(clpxfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_continuous_LASSOplus_XGBoost_validate"))
      }
    }else if (methods == "LASSO2plus_XGBoost"){
      cl2xfit = LASSO2plus_XGBtraining(data = tdat, biomks = varsCont, outcomeType = "continuous", Y = Ycont, outfile = paste0(outfileName,"_continuous_LASSO2plus_XGBoost"))
      if(!is.null(vdat)){
        pcl2xfit = XGBtraining_predict(cl2xfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_continuous_LASSO2plus_XGBoost_validate"))
      }
    }else{
      cl = LASSO2(data = tdat, biomks = varsCont, outcomeType = "continuous", Y = Ycont, outfile = paste0(outfileName,"_continuous_LASSO2"))
      clr = LASSO2_reg(data = tdat, biomks = varsCont, outcomeType = "continuous", Y = Ycont, outfile = paste0(outfileName,"_continuous_LASSO2reg"))
      cfit = LASSO_plus(data = tdat, biomks = varsCont, outcomeType = "continuous", Y = Ycont, outfile = paste0(outfileName,"_continuous_LASSOplus"))
      c2fit = LASSO2plus(data = tdat, biomks = varsCont, outcomeType = "continuous", Y = Ycont, outfile = paste0(outfileName,"_continuous_LASSO2plus"))
      cxfit = XGBtraining(data = tdat, biomks = varsCont, outcomeType = "continuous", Y = Ycont, outfile = paste0(outfileName,"_continuous_XGBoost"))
      clxfit = LASSO2_XGBtraining(data = tdat, biomks = varsCont, outcomeType = "continuous", Y = Ycont, outfile = paste0(outfileName,"_continuous_LASSO2_XGBoost"))
      clpxfit = LASSO_plus_XGBtraining(data = tdat, biomks = varsCont, outcomeType = "continuous", Y = Ycont, outfile = paste0(outfileName,"_continuous_LASSOplus_XGBoost"))
      cl2xfit =  LASSO2plus_XGBtraining(data = tdat, biomks = varsCont, outcomeType = "continuous", Y = Ycont, outfile = paste0(outfileName,"_continuous_LASSO2plus_XGBoost"))
      if(!is.null(vdat)){
        pcl = LASSO2_predict(cl, newdata = vdat, newY = TRUE, outfile =  paste0(outfileName,"_continuous_LASSO2_validate"))
        pclr = rms_model(clr$fit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_continuous_LASSO2reg_validate"))
        pcfit = rms_model(cfit$fit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_continuous_LASSOplus_validate"))
        pc2fit = rms_model(c2fit$fit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_continuous_LASSO2plus_validate"))
        pcxfit = XGBtraining_predict(cxfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_continuous_XGBoost_validate"))
        pclxfit = XGBtraining_predict(clxfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_continuous_LASSO2_XGBoost_validate"))
        pclpxfit = XGBtraining_predict(clpxfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_continuous_LASSOplus_XGBoost_validate"))
        pcl2xfit = XGBtraining_predict(cl2xfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_continuous_LASSO2plus_XGBoost_validate"))
      }
    }
  }

  
  ## time to event
  if(!is.null(time) & !is.null(event)){
    if (methods == "LASSO2"){
      tl = LASSO2(data = tdat, biomks = varsSurvival, outcomeType = "time-to-event", time = time, event = event, outfile = paste0(outfileName,"_survival_LASSO2"))
      if(!is.null(vdat)){
        ptl = LASSO2_predict(tl, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_survival_LASSO2_validate"))
      }
    }else if (methods == "LASSO2_reg"){
      tlr = LASSO2_reg(data = tdat, biomks = varsSurvival, outcomeType = "time-to-event", time = time, event = event, outfile = paste0(outfileName,"_survival_LASSO2reg"))
      if(!is.null(vdat)){
        ptlr = rms_model(tlr$fit, data = tdat, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_survival_LASSO2reg_validate"))
      }
    }else if (methods == "LASSO_plus"){
      tfit = LASSO_plus(data = tdat, biomks = varsSurvival, outcomeType = "time-to-event", time = time, event = event, outfile = paste0(outfileName,"_survival_LASSOplus"))
      if(!is.null(vdat)){
        ptfit = rms_model(tfit$fit, data = tdat, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_survival_LASSOplus_validate"))
      }
    }else if (methods == "LASSO2plus"){
      t2fit = LASSO2plus(data = tdat, biomks = varsSurvival, outcomeType = "time-to-event", time = time, event = event, outfile = paste0(outfileName,"_survival_LASSO2plus"))
      if(!is.null(vdat)){
        pt2fit = rms_model(t2fit$fit, data = tdat, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_survival_LASSO2plus_validate"))
      }
    }else if (methods == "XGBoost"){
      txfit = XGBtraining(data = tdat, biomks = varsSurvival, outcomeType = "time-to-event", time = time, event = event, outfile = paste0(outfileName,"_survival_XGBoost"))
      if(!is.null(vdat)){
        ptxfit = XGBtraining_predict(txfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_survival_XGBoost_validate"))
      }
    }else if (methods == "LASSO2_XGBoost"){
      tlxfit = LASSO2_XGBtraining(data = tdat, biomks = varsSurvival, outcomeType = "time-to-event", time = time, event = event, 
                                 outfile = paste0(outfileName,"_survival_LASSO2_XGBoost"))
      if(!is.null(vdat)){
        ptlxfit = XGBtraining_predict(tlxfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_survival_LASSO2_XGBoost_validate"))
      }
    }else if (methods == "LASSO_plus_XGBoost"){
      tlpxfit = LASSO_plus_XGBtraining(data = tdat, biomks = varsSurvival, outcomeType = "time-to-event", time = time, event = event, 
                                       outfile = paste0(outfileName,"_survival_LASSOplus_XGBoost"))
      if(!is.null(vdat)){
        ptlpxfit = XGBtraining_predict(tlpxfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_survival_LASSOplus_XGBoost_validate"))
      }
    }else if (methods == "LASSO2plus_XGBoost"){
      tl2xfit = LASSO2plus_XGBtraining(data = tdat, biomks = varsSurvival, outcomeType = "time-to-event", time = time, event = event, 
                                       outfile = paste0(outfileName,"_survival_LASSO2plus_XGBoost"))
      if(!is.null(vdat)){
        ptl2xfit = XGBtraining_predict(tl2xfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_survival_LASSO2plus_XGBoost_validate"))
      }
    }else{
      tl = LASSO2(data = tdat, biomks = varsSurvival, outcomeType = "time-to-event", time = time, event = event, outfile = paste0(outfileName,"_survival_LASSO2"))
      tlr = LASSO2_reg(data = tdat, biomks = varsSurvival, outcomeType = "time-to-event", time = time, event = event, outfile = paste0(outfileName,"_survival_LASSO2reg"))
      tfit = LASSO_plus(data = tdat, biomks = varsSurvival, outcomeType = "time-to-event", time = time, event = event, outfile = paste0(outfileName,"_survival_LASSOplus"))
      t2fit = LASSO2plus(data = tdat, biomks = varsSurvival, outcomeType = "time-to-event", time = time, event = event, outfile = paste0(outfileName,"_survival_LASSO2plus"))
      txfit = XGBtraining(data = tdat, biomks = varsSurvival, outcomeType = "time-to-event", time = time, event = event, outfile = paste0(outfileName,"_survival_XGBoost"))
      tlxfit = LASSO2_XGBtraining(data = tdat, biomks = varsSurvival, outcomeType = "time-to-event", time = time, event = event, 
                                 outfile = paste0(outfileName,"_survival_LASSO2_XGBoost"))
      tlpxfit = LASSO_plus_XGBtraining(data = tdat, biomks = varsSurvival, outcomeType = "time-to-event", time = time, event = event, 
                                       outfile = paste0(outfileName,"_survival_LASSOplus_XGBoost"))
      tl2xfit = LASSO2plus_XGBtraining(data = tdat, biomks = varsSurvival, outcomeType = "time-to-event", time = time, event = event, 
                                       outfile = paste0(outfileName,"_survival_LASSO2plus_XGBoost"))
      if(!is.null(vdat)){
        ptl = LASSO2_predict(tl, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_survival_LASSO2_validate"))
        ptlr = rms_model(tlr$fit, data = tdat, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_survival_LASSO2reg_validate"))
        ptfit = rms_model(tfit$fit, data = tdat, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_survival_LASSOplus_validate"))
        pt2fit = rms_model(t2fit$fit, data = tdat, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_survival_LASSO2plus_validate"))
        ptxfit = XGBtraining_predict(txfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_survival_XGBoost_validate"))
        ptlxfit = XGBtraining_predict(tlxfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_survival_LASSO2_XGBoost_validate"))
        ptlpxfit = XGBtraining_predict(tlpxfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_survival_LASSOplus_XGBoost_validate"))
        ptl2xfit = XGBtraining_predict(tl2xfit, newdata = vdat, newY = TRUE, outfile = paste0(outfileName,"_survival_LASSO2plus_XGBoost_validate"))
      }
    }
  }

  
  outs = list(bl,cl,tl,blr,clr,tlr,bfit, cfit,tfit, b2fit, c2fit, t2fit,bxfit, cxfit, txfit, blxfit, clxfit, tlxfit, blpxfit, clpxfit, tlpxfit,  bl2xfit, cl2xfit, tl2xfit,
              pbl,pcl,ptl,pblr,pclr,ptlr,pbfit, pcfit,ptfit, pb2fit, pc2fit, pt2fit,pbxfit, pcxfit, ptxfit, pblxfit, 
              pclxfit, ptlxfit, pblpxfit, pclpxfit, ptlpxfit,  pbl2xfit, pcl2xfit, ptl2xfit)
  names(outs) = c("LASSO2_binary", "LASSO2_cont", "LASSO2_surv",
                  "LASSO2reg_binary", "LASSO2reg_cont", "LASSO2reg_surv",
                  "LASSOplus_binary", "LASSOplus_cont", "LASSOplus_surv",
                  "LASSO2plus_binary", "LASSO2plus_cont", "LASSO2plus_surv",
                  "XGBoost_binary", "XGBoost_cont", "XGBoost_surv",
                  "LASSO2_XGBoost_binary", "LASSO2_XGBoost_cont", "LASSO2_XGBoost_surv",
                  "LASSOplus_XGBoost_binary", "LASSOplus_XGBoost_cont", "LASSOplus_XGBoost_surv",
                  "LASSO2plus_XGBoost_binary", "LASSO2plus_XGBoost_cont", "LASSO2plus_XGBoost_surv",
                  "pred_LASSO2_binary", "pred_LASSO2_cont", "pred_LASSO2_surv",
                  "pred_LASSO2reg_binary", "pred_LASSO2reg_cont", "pred_LASSO2reg_surv",
                  "pred_LASSOplus_binary", "pred_LASSOplus_cont", "pred_LASSOplus_surv",
                  "pred_LASSO2plus_binary", "pred_LASSO2plus_cont", "pred_LASSO2plus_surv",
                  "pred_XGBoost_binary", "pred_XGBoost_cont", "pred_XGBoost_surv",
                  "pred_LASSO2_XGBoost_binary", "pred_LASSO2_XGBoost_cont", "pred_LASSO2_XGBoost_surv",
                  "pred_LASSOplus_XGBoost_binary", "pred_LASSOplus_XGBoost_cont", "pred_LASSOplus_XGBoost_surv",
                  "pred_LASSO2plus_XGBoost_binary", "pred_LASSO2plus_XGBoost_cont", "pred_LASSO2plus_XGBoost_surv"
                  )
  outs = Filter(Negate(is.null), outs)
  return(outs)

}
  

