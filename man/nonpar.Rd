% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonpar.R
\name{nonpar}
\alias{nonpar}
\title{Nonparametric Random Graph Testing
Wrapper function to take in two graphs and run the nonpar
test with alignment; however it may be more useful to use nonpar.test
directly after alignment}
\usage{
nonpar(A, B, d = 5, nsims = 100, Q_init = NULL, lambda_init = 0.5,
  lambda_final = 0.01, eta = 0.5, eps = 0.01, numReps = 20,
  match = FALSE, alpha = 0.05)
}
\arguments{
\item{A}{An adjacency matrix}

\item{B}{another adjacency matrix}

\item{d}{the (known) dimension}

\item{nsims}{the number of permutations to run}

\item{Q_init}{initial guess for orthogonal matrix}

\item{lambda_init}{initial value of lambda for regularized sinkhorn problem}

\item{lambda_final}{final value of lambda}

\item{eta}{decay rate of the lambdas}

\item{eps}{tolerance}

\item{numReps}{reps in the iterative sinkhorn divergence problem}

\item{match}{If true, run a optimal-transport procrustes}

\item{alpha}{the type one error}
}
\description{
Nonparametric Random Graph Testing
Wrapper function to take in two graphs and run the nonpar
test with alignment; however it may be more useful to use nonpar.test
directly after alignment
}
