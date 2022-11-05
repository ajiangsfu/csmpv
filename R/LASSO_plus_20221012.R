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
#' 
#' wrote code first, later, update the above section

LASSO_plus= function(data = NULL, standardization = FALSE, columnWise = TRUE, biomks = NULL, allmks = FALSE, outputname = "test",
                     outcomeType = c("binary","continuous","time-to-event"), Y = NULL, time = NULL, event = NULL, outfile = "someName"){

  
  if(is.null(data)){
    stop("Please input a data set")
  }
  if(standardization){
    data = standardize(data, byrow = !columnWise)
  }
  
  aout = NA
  alls = NA
  
  if(outcomeType == "binary"){
    ## should call LASSO_plus for binary outcome, which I should already have
 
    
  }else if(outcomeType == "continuous"){
    ## I do not think that I have this function, but should be easily implemented
 
    
  }else if(outcomeType == "time-to-event"){
    
    ## should call LASSO_plus for binary outcome, which I should already have
    
  }else{
    stop("Please select the correct outcome type")
  }
  
  ## each output of a biomk is a list, which contains two items: fit object and fit coef
  
  ## combine all models together for forest plot
  fitout = lapply(aout, function(xx){
    xx= xx[[1]]
    xx = forestmodel::forest_model(xx, format_options = forest_model_format_options(text_size= 4, point_size = 4)) +
      theme(axis.text.x = element_text(size=4))
    return(xx)
  }) 
  
  kk = length(fitout)
  
  outplot = paste0(outfile, "3.pdf")
  pdf(outplot, height = kk, width = 7)
  ggarrange(plotlist = fitout, ncol=1, nrow=kk)
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
  aplot = paste("allMarks",outplot, sep="_")
  pdf(aplot, height = ceiling(kk/2), width = 7)
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
#outcomeType = "binary"
#Y = "CODE_postBMTFFS"

## test for cont Y
Y = "postBMTFFS"
outcomeType = "continuous"

## test for time to event
outcomeType = "time-to-event"
time = "postBMTFFS"
event = "CODE_postBMTFFS"
## all codes are working now! of course, I need to use the code file with latested date, if there are multiple files for same functions

### in the very end,before I package the file, I should start a new folder to save final version of code 
### and replace the whole folder if I want to make any changes






