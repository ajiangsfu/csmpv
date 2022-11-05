#' LASSO_plus for variable selection 
#' @description This is the wrap up function to select variables that are associated with an outcome variable based on my new algorithm: LASSO_plus. 
#' @details This function is to use LASSO_plus algorithm to select variables that are assoicated with on an outcome variable for a given data set.
#' An outcome variable could be a binary, continuous, or time-to-event variable. LASSO_plus selects variables based on LASSO, single variable regession, 
#' and stepwise regression.  
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

LASSO_plus = function(data = NULL, standardization = FALSE, columnWise = TRUE, biomks = NULL, outcomeType = c("binary","continuous","time-to-event"), 
                     Y = NULL, time = NULL, event = NULL, topN = 10, outfile = "someName", height = 6){

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
  aplot = paste("LASSO_plus_varaibleSelection",outfile, sep="_")
  pdf(aplot, height = height, width = 7)
  forestmodel::forest_model(xx, format_options = forest_model_format_options(text_size= 4, point_size = 4)) +
    theme(axis.text.x = element_text(size=4))
  dev.off()
  
  acoe = alls[[2]]
  acoeOut = gsub("pdf", "csv", aplot)
  write.csv(acoe, acoeOut)
  
}


############# testing code ###########
## read in data first:
library(rstudioapi)
current_working_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
# "/Users/aijiang/Desktop/AJworking/SepOct2022/varSelectPred_SepOct2022/varSelectPred_package/varSelectPred/R"
setwd(current_working_dir)
dat = read.csv("../../exampleData/proteinPerc_ASCT1_postBMTFFS_49ids_20220331.csv", header = T, row.names = 1, stringsAsFactors = F)
tmp = grep("Percent$", colnames(dat))
vars = colnames(dat)[tmp]

### test:



################################## stop here on 20221028#################
## got problem with this wrap-up function, fix it one next Monday: 20221031
## tried for binary outcome, and this is what I got:
# Error in forest_model_format_options(text_size = 4, point_size = 4) : 
#   could not find function "forest_model_format_options"
# In addition: There were 50 or more warnings (use warnings() to see the first 50)
###################################################



## binary
LASSO_plus(data = dat, biomks = vars, Y =  "CODE_postBMTFFS", outfile = "binaryOutcome_20221028")

## continuous
LASSO_plus(data = dat, biomks = vars, outcomeType = "continuous", Y =  "postBMTFFS", outfile = "continousOutcome_20221028")

## time to event
LASSO_plus(data = dat, biomks = vars, outcomeType = "time-to-event", time = "postBMTFFS",  event = "CODE_postBMTFFS", outfile = "timeToEvenOutcome_20221028")

