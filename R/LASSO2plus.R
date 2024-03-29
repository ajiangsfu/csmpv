#' @title Variable Selection and Modeling with LASSO2plus
#' 
#' @description
#' This function performs variable selection using the LASSO2plus algorithm and subsequently builds a model.
#' 
#' @details
#' The LASSO2plus algorithm begins with variable selection using LASSO2, typically involving multiple cross-validation-based LASSO regressions.
#' However, if only one or no variables are selected, the cross-validation results are ignored, and the algorithm ensures a minimum 
#' of two remaining variables through full-data lambda simulations. Additionally, it conducts variable selection through single-variable regression for each candidate variable.
#' The variables selected from both LASSO2 and single-variable approaches are then combined to perform traditional variable selection using stepwise regression. 

#' This function is designed to handle outcome variables of binary, continuous, or time-to-event type. Following variable selection, 
#' a model is constructed using standard R functions such as lm, glm, or coxph, depending on the type of outcome variable.
#' 
#' @param data A data matrix or a data frame, samples are in rows, and features/traits are in columns.
#' @param standardization A logic variable to indicate if standardization is needed before variable selection, the default is FALSE.
#' @param columnWise A logic variable to indicate if column wise or row wise normalization is needed, the default is TRUE, which is to do column-wise normalization. 
#'        This is only meaningful when "standardization" is TRUE.
#' @param biomks A vector of potential biomarkers for variable selection, they should be a subset of "data" column names. 
#' @param outcomeType Outcome variable type. There are three choices: "binary" (default), "continuous", and "time-to-event".  
#' @param Y Outcome variable name when the outcome type is either "binary" or "continuous". 
#' @param time Time variable name when outcome type is "time-to-event".
#' @param event Event variable name when outcome type is "time-to-event".
#' @param outfile A string for the output file including path if necessary but without file type extension. 
#' @param height An integer to indicate the forest plot height in inches
#' @importFrom utils write.csv
#' @author Aixiang Jiang
#' @return A list is returned: 
#' \item{fit}{A model with selected variables for the given outcome variable}
#' \item{outplot}{A forest plot}
#' @references 
#'   Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent (2010), Journal of Statistical Software, Vol. 33(1), 1-22, doi:10.18637/jss.v033.i01.
#'   
#'   Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol. 39(5), 1-13, doi:10.18637/jss.v039.i05.
#'   
#'   Hastie, T. J. and Pregibon, D. (1992) Generalized linear models. Chapter 6 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
#'   
#'   Therneau, T., Grambsch, P., Modeling Survival Data: Extending the Cox Model. Springer-Verlag, 2000.
#'   
#'   Kassambara A, Kosinski M, Biecek P (2021). survminer: Drawing Survival Curves using 'ggplot2'_. R package version 0.4.9,
#'         <https://CRAN.R-project.org/package=survminer>.
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

#' # The function can work with three different outcome types. 
#' # Here, we use continuous as an example:
#' c2fit = LASSO2plus(data = tdat, biomks = Xvars,
#'                    outcomeType = "continuous", Y = "Age",
#'                    outfile = paste0(temp_dir, "/continuousLASSO2plus"))
#' # You might save the files to the directory you want.
#' 
#' # To delete the "temp_dir", use the following:
#' unlink(temp_dir)

#' @export

LASSO2plus = function(data = NULL, standardization = FALSE, columnWise = TRUE, biomks = NULL, outcomeType = c("binary","continuous","time-to-event"), 
                      Y = NULL, time = NULL, event = NULL, outfile = "nameWithPath", height = 6){
  
  if(is.null(data)){
    stop("Please input a data set")
  }
  if(standardization){
    data = standardize(data, byrow = !columnWise)
  }
  
  alls = NA
  
  outcomeType = outcomeType[1]
  
  if(outcomeType == "binary"){
    alls = LASSO2plus_binary(data, biomks, Y, outfile = outfile)
    
  }else if(outcomeType == "continuous"){
    alls = LASSO2plus_continuous(data, biomks, Y, outfile = outfile)
    
  }else if(outcomeType == "time-to-event"){
    alls = LASSO2plus_timeToEvent (data, biomks,  time, event, outfile = outfile)
    
  }else{
    stop("Please select the correct outcome type")
  }
  
  xx= alls[[1]]
  aplot = paste0(outfile,"LASSO2plus_varaibleSelection.pdf")
  #pdf(aplot, height = height, width = 7)
  newplot = forestmodel::forest_model(xx, format_options = forestmodel::forest_model_format_options(text_size= 4, point_size = 4)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=4))
  ##dev.off()
  ggpubr::ggexport(newplot, filename = aplot, height = height, width = 7)
  
  acoe = alls[[2]]
  acoeOut = gsub("pdf", "csv", aplot)
  write.csv(acoe, acoeOut)
  fit = alls[[1]]
  outs = list(fit, newplot)
  names(outs) = c("fit", "outplot")
  return(outs) 
}
