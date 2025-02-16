% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gibbsCpp}
\alias{gibbsCpp}
\title{The gibbs function}
\usage{
gibbsCpp(N, a, b, n, burn)
}
\arguments{
\item{N}{size}

\item{a}{the param of beta}

\item{b}{the param of beta}

\item{n}{the size of x}

\item{burn}{start from burn.}
}
\value{
gibbs result
}
\description{
A function is used in homework
}
\examples{
\dontrun{
N <- 10000
a <- 2
b <- 3
n <- 10
burn <- 1000
result <- gibbsCpp(N, a, b, n, burn)
}
}
