% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_roc.R
\name{plot_roc}
\alias{plot_roc}
\title{Plot IPCW ROC curve}
\usage{
plot_roc(time, delta, marker, wts = NULL, tao, method)
}
\arguments{
\item{time}{vector of (censored) event-times}

\item{delta}{vector of event indicators at the corresponding value of the vector time. Censored observations must be denoted by the value 0.}

\item{marker}{vector of the marker values for which we want to compute the time-dependent ROC curve. the function assumes that larger values of the marker are associated with higher risks of events}

\item{wts}{IPC weights}

\item{tao}{evaluation time point of interest}

\item{method}{IPCW or discard}
}
\value{
a plot IPCW ROC curve
}
\description{
Plot IPCW ROC curve for a prediction
}
