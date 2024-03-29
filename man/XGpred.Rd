% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XGpred.R
\name{XGpred}
\alias{XGpred}
\title{XGpred: Building Risk Classification Predictive Models using Survival Data}
\usage{
XGpred(
  data = NULL,
  varsIn = NULL,
  selection = FALSE,
  vsMethod = c("LASSO2", "LASSO2plus", "LASSO_plus"),
  time = NULL,
  event = NULL,
  nrounds = 5,
  probcut = 0.8,
  nclass = c(2, 3),
  topN = 10,
  outfile = "nameWithPath"
)
}
\arguments{
\item{data}{A data matrix or a data frame where samples are in rows and features/traits are in columns.}

\item{varsIn}{A vector of variables used for the prediction model.}

\item{selection}{Logical. Default is FALSE. If TRUE, three variable selection methods can be chosen.}

\item{vsMethod}{When "selection" is set to TRUE, three variable selection methods can be chosen, with LASSO2 as the default method. The other two methods
are "LASSO2plus" and "LASSO_plus."}

\item{time}{Time variable name.}

\item{event}{Event variable name.}

\item{nrounds}{The maximum number of boosting iterations.}

\item{probcut}{Probability cutoff for risk group classification. Default is set to 0.8.}

\item{nclass}{Number of risk groups. By default, it is 2; any samples not classified into high-risk groups are classified into the low-risk group. 
When 3 is chosen, samples are classified into low, middle, and high-risk groups.}

\item{topN}{An integer indicating how many variables to select if LASSO_plus is chosen as the variable selection method.}

\item{outfile}{A string for the output file, including the path if necessary but without a file type extension.}
}
\value{
A list is returned with the following seven items:
\item{ranks}{Ranks from XGboost and Cox}
\item{twoEnds}{High and low risk group samples identified by mean ranks from XGBoost and Cox models}
\item{weights}{Weights for each variables used in the model}
\item{modelPars}{Mean and standard error of model scores for each risk group}
\item{nclass}{Number of risk groups}
\item{XGpred_score}{Model XGpred score}
\item{XGpred_prob}{Empirical Bayesian probability based on model XGpred score}
\item{XGpred_prob_class}{Risk group classification based on XGpred_prob for the given probability cutoff}
\item{probcut}{Probability cutoff for risk group classification}
}
\description{
The XGpred function is designed to generate an empirical Bayesian-based binary risk classification model with survival data based on our novel XGpred algorithm, 
combining XGBoost and traditional survival analysis.
}
\details{
If variable selection is needed, three variable selection methods are provided. Either the given variable or the selected variable list is used to build both an XGBoost model
and a traditional Cox model. Risk scores for each model are calculated and ranked, then averaged for each sample. 
The top 1/3 of samples are defined as the high-risk group, while the bottom 1/3 of samples are defined as the low-risk group. 
The binary risk classification model is built based on these two risk groups using either the given variable or the selected variable list. 
The model is a linear combination of these variables, with weights defined as t values derived from the single-variable linear model
of each variable on the two groups. Finally, the classification is based on empirical Bayesian probabilities.
}
\examples{
# Load in data sets:
data("datlist", package = "csmpv")
tdat = datlist$training

# The function saves files locally. You can define your own temporary directory. 
# If not, tempdir() can be used to get the system's temporary directory.
temp_dir = tempdir()
# As an example, let's define Xvars, which will be used later:
Xvars = c("highIPI", "B.Symptoms", "MYC.IHC", "BCL2.IHC", "CD10.IHC", "BCL6.IHC")
# For given time-to-event outcome and Xvars, we can build up a binary risk classification:
 xgobj = XGpred(data = tdat, varsIn = Xvars, 
                time = "FFP..Years.", event = "Code.FFP", 
                outfile = paste0(temp_dir, "/XGpred"))
# You might save the files to the directory you want.

# To delete the temp_dir, use the following:
unlink(temp_dir)
}
\references{
Tianqi Chen and Carlos Guestrin (2016), "XGBoost: A Scalable Tree Boosting System", 22nd SIGKDD Conference on Knowledge Discovery and Data Mining, 2016, https://arxiv.org/abs/1603.02754
 
 Aoki T, Jiang A, Xu A et al.,(2023) Spatially Resolved Tumor Microenvironment Predicts Treatment Outcomes in Relapsed/Refractory Hodgkin Lymphoma. J Clin Oncol. 2023 Dec 19:JCO2301115. doi: 10.1200/JCO.23.01115. Epub ahead of print. PMID: 38113419.
}
\author{
Aixiang Jiang
}
