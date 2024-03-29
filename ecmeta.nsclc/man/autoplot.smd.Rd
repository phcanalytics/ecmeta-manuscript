% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smd.R
\name{autoplot.smd}
\alias{autoplot.smd}
\alias{autoplot.smd_list}
\title{Plot standardized mean differences}
\usage{
\method{autoplot}{smd}(object, method = "iptw_att")

\method{autoplot}{smd_list}(object, method = "iptw_att", id = NULL)
}
\arguments{
\item{object}{An object of the appropriate class.}

\item{method}{A character vector of length one indicating the propensity score
method used to generate the weights.}

\item{id}{Character vector of the same length as \code{object} to name each element
in the list of \code{smd} objects. Passed to the \code{.id} argument in
\code{dplyr::bind_rows()}.}
}
\value{
A \code{\link[ggplot2:ggplot]{ggplot2::ggplot}} object.
}
\description{
Quickly plot standardized mean differences. Plots are faceted by group.
}
\details{
If \code{method} is \code{NULL}, then all methods are compared in a single plot.
However, in this case it is infeasible to display each covariate in the plot
so only the (logit of the) propensity scores are compared. If \code{object} is of
class \code{smd_list}, then SMDs are compared (via the \code{colour} aesthetic) across
all objects in the list.
}
