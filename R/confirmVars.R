#' Known biomarker confirmation 
#' @description This is the wrap up function to confirm/validate known biomarkers for a given data set
#' @details This function is to confirm/validate if a single variable or a group of variables has or have effect on an outcome variable for a given data set.
#' An outcome variable could be a binary, continuous, or time-to-event variable. 
#' @param data A data matrix or a data frame, samples are in columns, and features/traits are in rows.
#' @param standardization A logic variable to indicate if standardization is needed biomarker confirmation/validation, the default value is FALSE.
#' @param columnWise A logic variable to indicate if column wise or row wise normalization is needed, the default is TRUE. 
#'        This is only meaningful when "standardization" is TRUE.
#' @param biomks A vector of biomarker names to be confirmed/validated, they should be a subset of "data" column names. Each of these biomarkers will be validated separately.
#' @param allmks A logic variable to indicate if all "biomks" should be analyzed together as well. The default is FALSE.
#' @param outcomeType Outcome variable type. There are three choices: "binary" (default), "continuous", and "time-to-event".  
#' @param Y Outcome varialbe name when the outcome type is either "binary" or "continuous". 
#' @param time Time variable name when outcome type is "time-to-event".
#' @param event Event variable name when outcome type is "time-to-event".
#' @keywords biomarker confirmation, biomaker validation
#' @author Aixiang Jiang
#' @return A list with three items is returned: 
#' \item{name}{content}
#' \item{}{}

#' @references 
#' 
#' @export

confirmVars = function(data = NULL, standardization = FALSE, columnWise = TRUE, biomks = NULL, allmks = FALSE, outputname = "test",
                       outcomeType = c("binary","continuous","time-to-event"), Y = NULL, time = NULL, event = NULL){
  if(is.null(data)){
    stop("Please input a data set")
  }
  if(standardization){
   data = standardize(data, byrow = !columnWise)
  }

  if(outcomeType == "binary"){
    ## check one variable at a time
    aout = sapply(biomks, function(aX){
      res = glmBinary(data,Y,aX)
    })
    
    ## combine all models together for forest plot
    
    
    ## combine all numbers together to write out
    
    ## if allmks is TRUE, run a multivariable model
    
    ## also write out number results and generate forest plot
    
  }else if(outcomeType == "continuous"){
    
  }else if(outcomeType == "time-to-event"){
    
  }else{
    stop("Please select the correct outcome type")
  }

}


############# testing code ###########
## read in data first:

library(rstudioapi)
current_working_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
# "/Users/aijiang/Desktop/AJworking/SepOct2022/varSelectPred_SepOct2022/varSelectPred_package/varSelectPred/R"
setwd(current_working_dir)
dat = read.csv("../../exampleData/proteinPerc_ASCT1_postBMTFFS_49ids_20220331.csv", header = T, row.names = 1, stringsAsFactors = F)
# 
# > dim(dat)
# [1]  49 127

## 
data = dat
tmp = grep("Percent$", colnames(dat))
tmp = colnames(dat)[tmp]
bn = 5
biomks = sample(tmp, size = bn)
allmks = TRUE
outcomeType = "binary"
Y = "CODE_postBMTFFS"

