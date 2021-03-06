% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prediction_discard.R
\name{prediction_discard}
\alias{prediction_discard}
\title{Prediction discarding censored observations}
\usage{
prediction_discard(train.data, test.data, xnam, tao, ens.library,
  tuneparams = NULL)
}
\arguments{
\item{train.data}{a data.frame with at least the following variables: 
event-times (censored) in the first column,
event indicator in the second column and 
covariates/features that the user potentially want to use in building the preodiction model.
Censored observations must be denoted by the value 0. Main event of interest is denoted by 1.}

\item{test.data}{a data.frame with the same variables and names that the train.data}

\item{xnam}{vector with the names of the covariates to be included in the model}

\item{tao}{evaluation time point of interest}

\item{tuneparams}{a list of tune parameters for each machine learning procedure. Name them as gam_param, lasso_param, randomforest_param, svm_param, bart_param, knn_param, nn_param.
Default values are the same used for the simulation.}
}
\value{
a list with the predictions of each machine learning algorithm and the AUC of each of them
}
\description{
Ad-hoc technique that discards censored observations from the analysis
}
