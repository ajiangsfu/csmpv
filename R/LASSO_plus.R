#' LASSO_plus is for variable selection, and model buildup with the most commonly used R functions: lm, glm, and coxph for different outcome type.
#' @description This is the wrap up function to select variables that are associated with an outcome variable based on my new algorithm: LASSO_plus, 
#'        and build up a model afterwards.
#' @details This function is to use LASSO_plus algorithm to select variables that are assoicated with on an outcome variable for a given data set.
#' An outcome variable could be a binary, continuous, or time-to-event variable. LASSO_plus selects variables based on LASSO, single variable regession, 
#' and stepwise regression. After variable selection, a model is built-up as well. 
#' @param data A data matrix or a data frame, samples are in columns, and features/traits are in rows.
#' @param standardization A logic variable to indicate if standardization is needed before variable selection, the default is FALSE.
#' @param columnWise A logic variable to indicate if column wise or row wise normalization is needed, the default is TRUE, which us tod do column-wise normalization. 
#'        This is only meaningful when "standardization" is TRUE.
#' @param biomks A vector of potential biomarkers for variable selection, they should be a subset of "data" column names. 
#' @param outcomeType Outcome variable type. There are three choices: "binary" (default), "continuous", and "time-to-event".  
#' @param Y Outcome varialbe name when the outcome type is either "binary" or "continuous". 
#' @param time Time variable name when outcome type is "time-to-event".
#' @param event Event variable name when outcome type is "time-to-event".
#' @param topN An interger to indicate how many variables we intend to select 
#' @param outfile A string for the output file including path if necessary but without file type extension. 
#' @param height An integer to indicate the forest plot height in inches
#' @keywords variable selection
#' @author Aixiang Jiang
#' @return A model is returned: 
#' \item{fit}{a model with selected variables for the given outcome variable}
#' @references 
#' 
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
  aplot = paste0("LASSO_plus_varaibleSelection_",outfile,".pdf")
  #pdf(aplot, height = height, width = 7)
  forestmodel::forest_model(xx, format_options = forestmodel::forest_model_format_options(text_size= 4, point_size = 4)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=4))
  ##dev.off()
  ggplot2::ggsave(aplot, height = height, width = 7)
  
  acoe = alls[[2]]
  acoeOut = gsub("pdf", "csv", aplot)
  write.csv(acoe, acoeOut)
  fit = alls[[1]]
  return(fit) 
}

