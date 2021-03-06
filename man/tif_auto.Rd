% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estheta.R
\name{tif_auto}
\alias{tif_auto}
\title{Item Information Function for IRT 1~3PLM}
\usage{
tif_auto(
  para,
  min = -6,
  max = 6,
  length.out = 301,
  D = 1.702,
  labs = F,
  simplify = T
)
}
\arguments{
\item{para}{the item parameter data.frame estimated by \code{\link{estip}}}

\item{min}{the slope parameter}

\item{max}{the location parameter}

\item{length.out}{desired length of the sequence.}

\item{D}{a scale constant}

\item{labs}{logical. If TRUE, result data.frame has theta column.}

\item{simplify}{logical. If TRUE, result is vactor not data.frame.}
}
\description{
This function returns data.frame of iif in selected sequence.
}
