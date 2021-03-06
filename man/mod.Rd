% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{mod}
\alias{mod}
\title{give the modulus of a number}
\usage{
mod(value, modulo)
}
\arguments{
\item{value}{A number from which you want to compute the modulo}

\item{modulo}{The value of the modulus.}
}
\value{
The remainder of the division of value by modulo
}
\description{
This function return the modulus of a number. The difference with
the R built-in function is that the value of n\%n is n instead of 0.
}
\examples{
mod(15,12)
mod(12,12)
}
