#' @title LASSO_plus Variable Selection and Modeling
#' 
#' @description This function performs variable selection using the LASSO_plus algorithm and builds a model afterward.
#' 
#' @details The LASSO_plus algorithm combines LASSO, single variable regression, and stepwise regression
#' to select variables associated with an outcome variable in a given dataset. The outcome variable can be binary, 
#' continuous, or time-to-event. After variable selection, a model is built using common R functions such as lm, glm, 
#' and coxph, depending on the outcome type.
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
#' @param topN An integer indicating the desired number of variables to be selected.
#' @param outfile A string representing the output file, including the path if necessary, but without the file type extension
#' @param height An integer to indicate the forest plot height in inches
#' @importFrom utils write.csv
#' 
#' @author Aixiang Jiang
#' @return A list is returned: 
#' \item{fit}{A model with selected variables for the given outcome variable}
#' \item{outplot}{A forest plot}
#' @references 
#'  Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent (2010), Journal of Statistical Software, Vol. 33(1), 1-22, doi:10.18637/jss.v033.i01.
#'  
#'  Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol. 39(5), 1-13, doi:10.18637/jss.v039.i05.
#'  
#'  Hastie, T. J. and Pregibon, D. (1992) Generalized linear models. Chapter 6 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
#'  
#'  Therneau, T., Grambsch, P., Modeling Survival Data: Extending the Cox Model. Springer-Verlag, 2000.
#'  
#'  Kassambara A, Kosinski M, Biecek P (2021). survminer: Drawing Survival Curves using 'ggplot2', R package version 0.4.9,
#'         <https://CRAN.R-project.org/package=survminer>.
#'  
#'  Aoki T, Jiang A, Xu A et al.,(2023) Spatially Resolved Tumor Microenvironment Predicts Treatment Outcomes in Relapsed/Refractory Hodgkin Lymphoma. J Clin Oncol. 2023 Dec 19:JCO2301115. doi: 10.1200/JCO.23.01115. Epub ahead of print. PMID: 38113419.

#'       
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
#' # Here, we use binary as an example:
#' bfit = LASSO_plus(data = tdat, biomks = Xvars, Y = "DZsig", topN = 5,
#'                   outfile = paste0(temp_dir, "/binaryLASSO_plus"))
#' # You might save the files to the directory you want.
#' 
#' # To delete the "temp_dir", use the following:
#' unlink(temp_dir)

#' @export

LASSO_plus = function(data = NULL, standardization = FALSE, columnWise = TRUE, biomks = NULL, outcomeType = c("binary","continuous","time-to-event"), 
                      Y = NULL, time = NULL, event = NULL, topN = 10, outfile = "nameWithPath", height = 6){
  
  if(is.null(data)){
    stop("Please input a data set")
  }
  if(standardization){
    data = standardize(data, byrow = !columnWise)
  }
  
  alls = NA
  
  outcomeType = outcomeType[1]
  
  if(outcomeType == "binary"){
    alls = LASSO_plus_binary(data, biomks, Y, topN)
    
  }else if(outcomeType == "continuous"){
    alls = LASSO_plus_continuous(data, biomks, Y, topN)
    
  }else if(outcomeType == "time-to-event"){
    alls = LASSO_plus_timeToEvent (data, biomks,  time, event, topN)
    
  }else{
    stop("Please select the correct outcome type")
  }
  
  xx= alls[[1]]
  aplot = paste0(outfile,"_LASSO_plus_varaibleSelection.pdf")
  newplot = forestmodel::forest_model(xx, format_options = forestmodel::forest_model_format_options(text_size= 4, point_size = 4)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=4))
  ggpubr::ggexport(newplot, filename = aplot, height = height, width = 7)
  
  acoe = alls[[2]]
  acoeOut = gsub("pdf", "csv", aplot)
  write.csv(acoe, acoeOut)
  fit = alls[[1]]
  outs = list(fit, newplot)
  names(outs) = c("fit", "outplot")
  return(outs) 
}
