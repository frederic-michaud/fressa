#' Plot the evolution of the frequency of haplotypes
#'
#' given a matrix of frequency returned by the function `compute.frequency.evolution`
#' and the associated genome, plot the frequency of all possible haplotype through time
#'
#' @param genome A S4 object of type genome
#' @param freqs a matrix of frequency as returned by the function `compute.frequency.evolution`
#' @examples
#' locus1 = create.locus(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1))
#' locus2 = create.locus(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#' freqs <- compute.frequency.evolution(genome)
#' plot.haplotype.frequency(genome, freqs)
#' @export

plot.haplotype.frequency <- function(genome,freqs){
  haplotype.frequency <- get.haplotype.frequency(genome,freqs)
  max.freq <- max(haplotype.frequency)
  palette = get.palette(get.nb.haplotype(genome))
  plot(haplotype.frequency[1,],type="l",ylim=c(0,1.2*max.freq),col=palette[1],xlab = "Generation",ylab="frequency")
  for(haplotype in 2:get.nb.haplotype(genome)){
    lines(haplotype.frequency[haplotype,],col=palette[haplotype])
  }
  legend("topright",legend=get.haplotype.names(genome),lty = rep(1,get.nb.haplotype(genome)),col=palette)
}

#' Plot the evolution of the frequency of genotype
#'
#' given a matrix of frequency returned by the function `compute.frequency.evolution`
#' and the associated genome, plot the frequency of all possible genotype through time
#'
#' @param genome A S4 object of the type genome
#' @param freqs a matrix of frequency as returned by the function `compute.frequency.evolution`
#' @examples
#' locus1 = create.locus(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1))
#' locus2 = create.locus(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#' freqs <- compute.frequency.evolution(genome)
#' plot.genotype.frequency(genome, freqs)
#' @export

plot.genotype.frequency <- function(genome,freqs){
  max.freq <- max(freqs)
  palette=get.palette(get.nb.genotype(genome))
  plot(freqs[1,],type="l",ylim=c(0,1.2*max.freq),col=palette[1],xlab = "Generation",ylab="frequency")
  for(genotype in 2:get.nb.genotype(genome)){
    lines(freqs[genotype,],col=palette[genotype])
  }
  legend("topright",legend=get.genotype.names(genome),lty = rep(1,get.nb.genotype(genome)),col=palette)
}

#' Plot the evolution of the frequency of an allele
#'
#' given a matrix of frequency returned by the function `compute.frequency.evolution`
#' and the associated genome, plot the frequency of all possible allele at a given locus.
#'
#' @param genome A S4 object of type genome
#' @param freqs a matrix of frequency as returned by the function `compute.frequency.evolution`
#' @param locus.position the index of the locus from which we want to plot the allele frequency
#' @examples
#' locus1 = create.locus(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1))
#' locus2 = create.locus(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#' freqs <- compute.frequency.evolution(genome)
#' plot.allele.frequency(genome, freqs,1)
#' plot.allele.frequency(genome, freqs,2)
#' @export

plot.allele.frequency <- function(genome,freqs,locus.position){
  allele.frequency <- get.allele.frequency(genome,freqs,locus.position)
  max.freq <- max(allele.frequency)
  allele.number <- get.nb.alleles.per.locus(genome)[locus.position]
  palette <- get.palette(allele.number)
  plot(allele.frequency[1,],type="l",ylim=c(0,1.2*max.freq),col=palette[1],xlab = "Generation",ylab="frequency")
  for(allele in 2:allele.number){
    lines(allele.frequency[allele,],col=palette[allele])
  }
  legend("topright",legend=get.allele.name(genome,locus.position),lty = rep(1,allele.number),col=palette)
}

#' Plot the marginal fitness of all haplotype in the population
#'
#' given a matrix of frequency returned by the function `compute.frequency.evolution`
#' and the associated genome, plot the evolution of marginal fitness. The marginal fitness i
#' s defined as the mean fitness of
#' individual carrying this haplotype weighted by the frequency of those individuals.
#'
#' @param genome A S4 object of type genome
#' @param freqs a matrix of frequency as returned by the function `compute.frequency.evolution`
#' @examples
#' locus1 = create.locus(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1))
#' locus2 = create.locus(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#' freqs <- compute.frequency.evolution(genome)
#' plot.haplotype.marginal.fitness(genome, freqs)
#' @export


plot.haplotype.marginal.fitness <- function(genome,freqs){
  haplotype.marginal.fitness <- get.marginal.haplotype.fitness(genome,freqs)
  max.fit <- max(haplotype.marginal.fitness,na.rm = T)
  min.fit <- min(haplotype.marginal.fitness,na.rm = T)
  haplotype.number <- get.nb.haplotype(genome)
  palette <- get.palette(haplotype.number)
  plot(haplotype.marginal.fitness[1,],type="l",ylim=c(min.fit/1.2,1.2*max.fit),col=palette[1],xlab = "Generation",ylab="fitness")
  for(haplotype in 2:haplotype.number){
    lines(haplotype.marginal.fitness[haplotype,],col=palette[haplotype])
  }
  legend("topright",legend=get.haplotype.names(genome),lty = rep(1,haplotype.number),col=palette)
}

#' Plot the marginal fitness of all allele from one locus
#'
#' given a matrix of frequency returned by the function `compute.frequency.evolution`
#' and the associated genome, plot the evolution of the marginal fitness of all allele
#' present in the population. The marginal fitness is defined as the mean fitness of
#' individual carrying this allele weighted by the frequency of those individuals.
#'
#' @param genome A S4 object of type genome
#' @param freqs a matrix of frequency as returned by the function `compute.frequency.evolution`
#' @param locus.position the index of the locus from which we want to plot the allele frequency
#' @examples
#' locus1 = create.locus(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1))
#' locus2 = create.locus(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#' freqs <- compute.frequency.evolution(genome)
#' plot.haplotype.marginal.fitness(genome, freqs)
#' @export


plot.allele.marginal.fitness <- function(genome,freqs,locus.position){
  allele.marginal.fitness <- get.marginal.allele.fitness(genome,freqs,locus.position)
  max.fit <- max(allele.marginal.fitness)
  min.fit <- min(allele.marginal.fitness)
  allele.number <- get.nb.alleles.per.locus(genome)[locus.position]
  palette <- get.palette(allele.number)
  plot(allele.marginal.fitness[1,],type="l",ylim=c(min.fit/1.2,1.2*max.fit),col=palette[1],xlab = "Generation",ylab="fitness")
  for(allele in 2:allele.number){
    lines(allele.marginal.fitness[allele,],col=palette[allele])
  }
  legend("topright",legend=get.allele.name(genome,locus.position),lty = rep(1,allele.number),col=palette)
}
