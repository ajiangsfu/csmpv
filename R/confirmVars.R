#' Biomarker confirmation 
#' @description This is the wrap up function to confirm/validate known biomarkers in a given data set
#' @details This function is to confirm/validate if a single variable or a group of variables has or have effect on an outcome variable for a given data set.
#' An outcome variable could be a binary, continuous, or a time-to-event variable. 
#' @param data A data matrix or a data frame, samples are in columns, and features/traits including outcome and biomarkers are in rows.
#' @param standardization A logic parameter to indicate if standardization is needed before biomarker confirmation/validation, the default value is FALSE.
#' @param columnWise A logic variable to indicate if column wise or row wise normalization is needed for standardization process, the default is TRUE. 
#'        This is only meaningful when "standardization" is TRUE.
#' @param biomks A vector of biomarker names to be confirmed/validated, they should be a subset of "data" column names. Each of these biomarkers will be 
#'        confirmed/validated separately, and all together if allmks is TRUE.
#' @param allmks A logic variable to indicate if all "biomks" should be analyzed together in addition to single variable analysis. The default is FALSE.
#' @param outcomeType Outcome variable type. There are three choices: "binary" (default), "continuous", and "time-to-event".  
#' @param Y Outcome varialbe name when the outcome type is either "binary" or "continuous". 
#' @param time Time variable name when outcome type is "time-to-event".
#' @param event Event variable name when outcome type is "time-to-event".
#' @param outfile A string for output files including path if necessary but without file type extension. 
#' @keywords biomarker confirmation, biomaker validation
#' @author Aixiang Jiang
#' 
#' @references 
#' @export

confirmVars = function(data = NULL, standardization = FALSE, columnWise = TRUE, biomks = NULL, allmks = FALSE, 
                       outcomeType = c("binary","continuous","time-to-event"), Y = NULL, time = NULL, event = NULL, outfile = "nameWithPath"){
  if(is.null(data)){
    stop("Please input a data set")
  }
  if(standardization){
   data = standardize(data, byrow = !columnWise)
  }

  aout = NA
  alls = NA
  
  outcomeType = outcomeType[1]
  
  if(outcomeType == "binary"){
    ## check one variable at a time
    aout = lapply(biomks, function(aX){
      res = glmBinary(data,Y,aX)
    })
    
    ## if allmks is TRUE, run a multivariable model
    if(allmks == TRUE){
      alls = glmBinary(data,Y,biomks)
    }
    
  }else if(outcomeType == "continuous"){
    ## check one variable at a time
    aout = lapply(biomks, function(aX){
      res = lmContinuous(data,Y,aX)
    })
    
    ## if allmks is TRUE, run a multivariable model
    if(allmks == TRUE){
      alls = lmContinuous(data,Y,biomks)
    }
    
  }else if(outcomeType == "time-to-event"){
    aout = lapply(biomks, function(aX){
      res = coxTimeToEvent(data, time, event, aX)
    })
    
    ## if allmks is TRUE, run a multivariable model
    if(allmks == TRUE){
      alls = coxTimeToEvent(data, time, event, biomks)
    }
    
  }else{
    stop("Please select the correct outcome type")
  }
  
  ## each output of a biomk is a list, which contains two items: fit object and fit coef
  
  ## combine all models together for forest plot
  fitout = lapply(aout, function(xx){
    xx= xx[[1]]
    xx = forestmodel::forest_model(xx, format_options = forestmodel::forest_model_format_options(text_size= 4, point_size = 4)) +
      theme(axis.text.x = element_text(size=4))
    return(xx)
  }) 
  
  kk = length(fitout)
  
  outplot = paste0(outfile, "3.pdf")
  pdf(outplot, height = kk, width = 7)
  ggpubr::ggarrange(plotlist = fitout, ncol=1, nrow=kk)
  dev.off()
  
  ## combine all numbers together to write out
  coeout = lapply(aout, function(xx){
    xx= xx[[2]]
    return(xx)
  }) 
  
  coes = do.call(rbind, coeout)
  ## should I keep Intercept for all models, maybe it is a good idea to keep to avoid confusion?
  
  if(rownames(coes)[1] == "1"){
    rownames(coes) = biomks
  }
  
  coeout = paste0(outfile, "3.csv")
  write.csv(coes, coeout)

  xx= alls[[1]]
  if(allmks){
    aplot = paste("allMarks",outplot, sep="_")
    pdf(aplot, height = ceiling(kk/2), width = 7)
    forestmodel::forest_model(xx, format_options = forestmodel::forest_model_format_options(text_size= 4, point_size = 4)) +
      theme(axis.text.x = element_text(size=4))
    dev.off()
    acoe = alls[[2]]
    acoeOut = gsub("pdf", "csv", aplot)
    write.csv(acoe, acoeOut)
  }

}

