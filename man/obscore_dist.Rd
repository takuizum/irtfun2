% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mext-contract.R
\name{obscore_dist}
\alias{obscore_dist}
\title{A function calculates IRT observed score using recursion formula.}
\usage{
obscore_dist(
  theta,
  a,
  b,
  c,
  D = 1.702,
  name = "test",
  color = "cyan",
  output = 1
)
}
\arguments{
\item{theta}{a vector of theta estimator EAP, MAP or MLE...}

\item{a}{a slope parameter.}

\item{b}{a location parameter}

\item{c}{a lower asymptote parameter}

\item{D}{a factor constant}

\item{name}{a plot title}

\item{color}{a plot color.}

\item{output}{int. if 1 score vector, if 2 cumulative distribution plot.}
}
\description{
A function calculates IRT observed score using recursion formula.
}
