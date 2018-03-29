#' Compute the evolution of the frequency of various genotype
#'
#' Given a genome, this function simulate the evolution of the frequency
#' of all the possible genotype that exist
#' @return A matrix containing the frequencies of each genotype at each generation
#' @param genome a list of loci
#' @param initial.frequency The initial frequency of the various genotype. If NULL is given
#' the initial frequencies will all be set to the same value
#' @param generations The number of generation that have to be computed (including the first one)
#' @examples
#' locus1 = data.frame(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1),fitness.male=c(1,1),fitness.female=c(1,1))
#' locus2 = data.frame(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = list(locus1,locus2)
#' freqs <- compute.frequency.evolution(genome)
#' plot(freqs[,1],type="l",col=1,ylim=c(0,.5),xlim = c(0,100))
#' lines(freqs[,2],type="l",col=2)
#' lines(freqs[,3],type="l",col=3)
#'
#'


compute.frequency.evolution <- function(genome,initial.frequency = NULL,generations = 25)
{
  nb.genotype <- get.nb.genotype(genome)
  if(is.null(initial.frequency)){
    initial.frequency <- rep(1/nb.genotype,nb.genotype)
  }
  if(length(initial.frequency)!=nb.genotype) stop("The given initial frequency is not the correct size")

  freqs <- rbind(initial.frequency,simulate.frequency(genome,initial.frequency)) #we compute the first generation here to be able to use tail
  for(generation in 1:(generations - 2))
  {
    freqs <- rbind(freqs,simulate.frequency(genome,tail(freqs,1))) #and the next 1000 generation
  }
  return(freqs)
}
