% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zero_gen.R
\name{generate_one_each_a_b_k}
\alias{generate_one_each_a_b_k}
\title{Generates one sample for each triplet of (a, b) with a fixed k under the abk parametrization.}
\usage{
generate_one_each_a_b_k(A, B, k)
}
\arguments{
\item{A}{A vector of parameter \code{a}s.}

\item{B}{A vector of parameter \code{b}s, of equal length as \code{A}.}

\item{k}{A positive number, the \code{k} parameter.}
}
\value{
A vector of the same size as \code{log_odds}.
}
\description{
Generates one sample for each pair of (a,b) with a fixed k under the abk parametrization.
}
\details{
Draws one example from the 1-d Hurdle with parameters \code{A[i]}, \code{B[i]}, \code{k} for each \code{i}.
}
\examples{
generate_one_each_a_b_k(-3:3, 0:6/10, 1)
}
