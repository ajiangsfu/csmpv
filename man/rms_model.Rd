% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rms_model.R
\name{rms_model}
\alias{rms_model}
\title{rms_model is to use rms package to build up a predictive model}
\usage{
rms_model(afit, data = NULL, newdata = NULL, outfile)
}
\arguments{
\item{afit}{A LASSO_plus object returned by LASSO_plus, or a model fit from any other R functions/programs}

\item{data}{A data frame that used to get afit, only this is only needed when the outcome variable is time to event, otherwise data is extracted from afit.}

\item{newdata}{A new data frame for prediction}

\item{outfile}{A string to indicate output file name without file type but should include all path information}
}
\description{
This is the function to build up a predictive model based on rms package after variable selection with LASSO_plus or a model from other R programs/fucntions.
}
\details{
This function is to build up a predictive model after variable selection with LASSO_plus or a model from other R programs/fucntions. 
      R package rms allows us to easily generate calibration (bootstrap) and nomograph plots, calculate C-index when necessary, predict outcome for a new data set. 
      If a rms based model is not needed, the model from LASSO_plus object can be directly used as well.
      For both of these models, a new data set should be comparable to the data set used for variable selection.
      If a new data set is not given, the training data set is used for prediction. In this function, validation and calibration are based on bootstrap.
}
\references{

}
\author{
Aixiang Jiang
}
\keyword{model}
\keyword{predictive}
