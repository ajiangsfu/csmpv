#' @title Biomarker Confirmation Function
#' @description This function confirms and validates known biomarkers in a given dataset.
#' @details Use this function to assess whether individual variables or groups of variables
#' have an impact on an outcome variable within a dataset. The outcome variable can be binary,
#' continuous, or time-to-event. Note that this function is not intended for model confirmation,
#' as it doesn't incorporate coefficients from previous research.
#' @import ggplot2
#' @import ggpubr
#' @importFrom utils write.csv
#' 
#' @param data A data matrix or data frame with samples in rows and features/traits (including outcome and biomarkers) in columns.
#' @param standardization Logical; indicates if standardization is needed before biomarker confirmation/validation. Default is FALSE.
#' @param columnWise Logical; indicates if column-wise or row-wise normalization is needed for standardization. Default is TRUE.
#' @param biomks A vector of biomarker names to confirm/validate. Subset of column names in the data parameter.
#' @param outcomeType The type of the outcome variable. It has three choices: "binary" (default), "continuous", and "time-to-event".
#' @param Y The outcome variable name when the outcome type is either "binary" or "continuous".
#' @param time The time variable name when the outcome type is "time-to-event".
#' @param event The event variable name when the outcome type is "time-to-event".
#' @param outfile A string representing the output file, including the path if necessary, but without the file type extension
#'
#' @return A list containing:
#' \item{fit}{A model with selected variables for the given outcome variable.}
#' \item{allplot}{A list of plots generated during the confirmation/validation process.}
#' There might be extra plots in the list for time-to-event outcome
#' @author Aixiang Jiang
#' @references 
#'  Hastie, T. J. and Pregibon, D. (1992) Generalized linear models. Chapter 6 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
#'  
#'  Therneau, T., Grambsch, P., Modeling Survival Data: Extending the Cox Model. Springer-Verlag, 2000.
#'  
#'  Kassambara A, Kosinski M, Biecek P (2021). survminer: Drawing Survival Curves using 'ggplot2', R package version 0.4.9,
#'         <https://CRAN.R-project.org/package=survminer>.

#' @examples
#' # Load in data sets:
#' data("datlist", package = "csmpv")
#' tdat = datlist$training
#' 
#' # The confirmVars function saves files locally. You can define your own temporary directory. 
#' # If not, tempdir() can be used to get the system's temporary directory.
#' temp_dir = tempdir()
#' 
#' # As an example, let's define Xvars, which will be used later:
#' Xvars = c("highIPI","B.Symptoms", "MYC.IHC", "BCL2.IHC", "CD10.IHC", "BCL6.IHC")
#' 
#' # confirmVars can work with three different outcome types. 
#' # Here, we use binary as an example:
#' bconfirm = confirmVars(data = tdat, biomks = Xvars, Y = "DZsig",
#'                         outfile = paste0(temp_dir, "/confirmBinary"))
#' # You might save the files to the directory you want.
#' 
#' # To delete the temp_dir, use the following:
#' unlink(temp_dir)

#' @export

confirmVars = function(data = NULL, standardization = FALSE, columnWise = TRUE, biomks = NULL,  
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
    
    alls = glmBinary(data,Y,biomks)
    
  }else if(outcomeType == "continuous"){
    ## check one variable at a time
    aout = lapply(biomks, function(aX){
      res = lmContinuous(data,Y,aX)
    })
    
    alls = lmContinuous(data,Y,biomks)

  }else if(outcomeType == "time-to-event"){
    aout = lapply(biomks, function(aX){
      res = coxTimeToEvent(data, time, event, aX, outfile)
    })
  
    alls = coxTimeToEvent(data, time, event, biomks, outfile)
    
  }else{
    stop("Please select the correct outcome type")
  }
  
  ## each output of a biomk is a list, which contains two items: fit object and fit coef
  
  ## combine all models together for forest plot
  fitout = lapply(aout, function(xx){
    xx= xx[[1]]
    xx = forestmodel::forest_model(xx, 
                                   format_options = forestmodel::forest_model_format_options(text_size= 4, point_size = 4)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size=4))
    invisible(xx)
  }) 
  
  kk = length(fitout)
  
  outplot = paste0(outfile, ".pdf")
  gplot = ggpubr::ggarrange(plotlist = fitout, ncol=1, nrow=kk)
  ggpubr::ggexport(gplot, filename = outplot, height = kk+2, width = 7) 
  allplot = list(gplot)
  
  ## combine all numbers together to write out
  coeout = lapply(aout, function(xx){
    xx= xx[[2]]
    return(xx)
  }) 
  
  coes = do.call(rbind, coeout)

  for(i in 1:length(biomks)){
    if(rownames(coes)[i] == as.character(i)){
      rownames(coes)[i] = biomks[i]
    }
  }

  coeout = paste0(outfile, ".csv")
  write.csv(coes, coeout)
  
  xx= alls[[1]]
  aplot = gsub("\\.pdf", "allMarks.pdf", outplot)
  gplot = forestmodel::forest_model(xx, 
                                    format_options = forestmodel::forest_model_format_options(text_size= 4, point_size = 4)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=4))
  ggpubr::ggexport(gplot, filename = aplot, height = ceiling(kk/2), width = 7)
  allplot = c(allplot, list(gplot))
  acoe = alls[[2]]
  acoeOut = gsub("pdf", "csv", aplot)
  write.csv(acoe, acoeOut)

  outs = list(xx, allplot)
  names(outs) = c("fit", "allplot")
  
  if(length(alls)>2){
    outs = c(outs, alls[[3]])
    names(outs)[3] = "extraPlot"
  }
  
  return (outs)

}
