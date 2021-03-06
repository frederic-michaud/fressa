% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.gamete.marginal.fitness}
\alias{plot.gamete.marginal.fitness}
\title{Plot the marginal fitness of all gamete in the population}
\usage{
plot.gamete.marginal.fitness(genome, freqs)
}
\arguments{
\item{genome}{A S4 object of type genome}

\item{freqs}{a matrix of frequency as returned by the function `compute.frequency.evolution`}
}
\description{
given a matrix of frequency returned by the function `compute.frequency.evolution`
and the associated genome, plot the evolution of marginal fitness. The marginal fitness i
s defined as the mean fitness of
individual carrying this gamete weighted by the frequency of those individuals.
}
\examples{
locus1 = create.locus(allele1=c(1,1),allele2 = c(1,2),sd = c(0,1))
locus2 = create.locus(allele1=  c(1,1,2),allele2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
genome = create.genome(locus=list(locus1,locus2))
freqs <- compute.frequency.evolution(genome)
plot.gamete.marginal.fitness(genome, freqs)
}
