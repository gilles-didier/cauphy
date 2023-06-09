% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimization.R
\name{initLambdaParameter}
\alias{initLambdaParameter}
\title{Initialization of the lambda parameter.}
\usage{
initLambdaParameter(tree, trait, disp_hat, tol = 0.1)
}
\arguments{
\item{trait}{named vector of traits at the tips.}

\item{disp_hat}{the previously estimated dispersion}

\item{tol}{the percentage of tip pairs to keep for the computation. Default to 0.1.}
}
\value{
The initial lambda parameter.
}
\description{
Initialize the lambda parameter.
}
\details{
Use linear combinations of cherries to get a first estimate of lambda.
}
\keyword{internal}
