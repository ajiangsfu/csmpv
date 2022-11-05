#' LASSO_plus for variable selection 
#' @description This is the wrap up function to select variables that are associated with an outcome variable based on my new algorithm: LASSO_plus. 
#' @details This function is to use LASSO_plus algorithm to select variables that are assoicated with on an outcome variable for a given data set.
#' An outcome variable could be a binary, continuous, or time-to-event variable. LASSO_plus selects variables based on LASSO, single variable regession, 
#' and stepwise regression.  
#' @param data A data matrix or a data frame, samples are in columns, and features/traits are in rows.
#' @param standardization A logic variable to indicate if standardization is needed biomarker confirmation/validation, the default value is FALSE.
#' @param columnWise A logic variable to indicate if column wise or row wise normalization is needed, the default is TRUE. 
#'        This is only meaningful when "standardization" is TRUE.
#' @param biomks A vector of potential biomarkers, they should be a subset of "data" column names. 
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

LASSO_plus= function(data, standardization, columnWise, biomks, allmks, outcomeType, Y, time, event){
  
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






