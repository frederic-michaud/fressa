#' Compute the evolution of the frequency of various genotype
#'
#' Given a genome, this function simulate the evolution of the frequency
#' of all the possible genotype that exist
#' @return A matrix containing the frequencies of each genotype at each generation
#' @param genome A S4 object of the type genome
#' @param initial.frequency The initial frequency of the various genotype. If NULL is given
#' the initial frequencies will all be set to the same value
#' @param generations The number of generation that have to be computed (including the first one)
#' @examples
#' locus1 = data.frame(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1),fitness.male=c(1,1),fitness.female=c(1,1))
#' locus2 = data.frame(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = list(locus1,locus2)
#' freqs <- compute.frequency.evolution(genome)
#'
#'


compute.frequency.evolution <- function(genome,initial.frequency = NULL,generations = 25)
{
  nb.genotype <- get.nb.genotype(genome)
  if(is.null(initial.frequency)){
    initial.frequency <- rep(1/nb.genotype,nb.genotype)
  }
  if(length(initial.frequency)!=nb.genotype) stop("The given initial frequency is not the correct size")

  freqs <- cbind(initial.frequency,simulate.frequency(genome,initial.frequency)) #we compute the first generation here to be able to use tail
  for(generation in 1:(generations - 2))
  {
    freqs <- cbind(freqs,simulate.frequency(genome,freqs[,ncol(freqs)])) #and the next 1000 generation
  }
  return(freqs)
}

#' Plot the evolution of the frequency of haplotypes
#'
#' given a matrix of frequency returned by the function `compute.frequency.evolution`
#' and the associated genome, plot the frequency of all possible haplotype through time
#'
#' @param genome A S4 object of the type genome
#' @param freqs a matrix of frequency as returned by the function `compute.frequency.evolution`
#' @examples
#' locus1 = create.locus(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1),fitness.male=c(1,1),fitness.female=c(1,1))
#' locus2 = create.locus(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#' freqs <- compute.frequency.evolution(genome)
#' plot.haplotype.frequency(genome, freqs)
#'

plot.haplotype.frequency <- function(genome,freqs){
  haplotype.frequency <- get.haplotype.frequency(genome,freqs)
  max.freq <- max(haplotype.frequency)
  plot(haplotype.frequency[1,],type="l",ylim=c(0,1.2*max.freq),col=1,xlab = "Generation",ylab="frequency")
  for(haplotype in 2:get.nb.haplotype(genome)){
    lines(haplotype.frequency[haplotype,],col=haplotype)
  }
  legend("topright",legend=get.haplotype.names(genome),lty = rep(1,get.nb.haplotype(genome)),col=1:get.nb.haplotype(genome))
}

get.haplotype.frequency <- function(genome,freqs)
{
  nb.generation <- ncol(freqs)
  haplotype.frequency <- matrix(0,ncol = nb.generation,nrow = get.nb.haplotype(genome))
  for (generation in 1:nb.generation){
    haplotype.frequency[,generation] <- get.haplotype.frequency.single.generation(genome,freqs[,generation])
  }
  row.names(haplotype.frequency) <- get.haplotype.names(genome)
  return(haplotype.frequency)
}

get.haplotype.frequency.single.generation <- function(genome,freqs)
{
  nb.haplotype <- get.nb.haplotype(genome)
  nb.genotype <- get.nb.genotype(genome)
  all.genotype <- genome@all.genotype
  #frequency is a matrix with the frequency of the various genome
  #row index indicate first haplotype and column index indicate second haplotype
  frequency <- matrix(0,nrow = nb.haplotype,ncol=nb.haplotype)
  for (genotype in 1:nb.genotype){
    frequency[all.genotype[genotype,1],all.genotype[genotype,2]] <- freqs[genotype]
  }
  sum.column <- colSums(frequency)/2
  sum.row <- rowSums(frequency)/2
  return(sum.column + sum.row)
}

#' Plot the evolution of the frequency of genotype
#'
#' given a matrix of frequency returned by the function `compute.frequency.evolution`
#' and the associated genome, plot the frequency of all possible genotype through time
#'
#' @param genome A S4 object of the type genome
#' @param freqs a matrix of frequency as returned by the function `compute.frequency.evolution`
#' @examples
#' locus1 = create.locus(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1),fitness.male=c(1,1),fitness.female=c(1,1))
#' locus2 = create.locus(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#' freqs <- compute.frequency.evolution(genome)
#' plot.genotype.frequency(genome, freqs)
#'


plot.genotype.frequency <- function(genome,freqs){
  max.freq <- max(freqs)
  plot(freqs[1,],type="l",ylim=c(0,1.2*max.freq),col=1,xlab = "Generation",ylab="frequency")
  for(genotype in 2:get.nb.genotype(genome)){
    lines(freqs[genotype,],col=genotype)
  }
  legend("topright",legend=get.genotype.names(genome),lty = rep(1,get.nb.genotype(genome)),col=1:get.nb.genotype(genome))
}

#' Plot the evolution of the frequency of an allele
#'
#' given a matrix of frequency returned by the function `compute.frequency.evolution`
#' and the associated genome, plot the frequency of all possible allele at a given locus.
#'
#' @param genome A S4 object of the type genome
#' @param freqs a matrix of frequency as returned by the function `compute.frequency.evolution`
#' @param locus.position the index of the locus from which we want to plot the allele frequency
#' @examples
#' locus1 = create.locus(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1),fitness.male=c(1,1),fitness.female=c(1,1))
#' locus2 = create.locus(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#' freqs <- compute.frequency.evolution(genome)
#' plot.allele.frequency(genome, freqs,1)
#' plot.allele.frequency(genome, freqs,2)
#'

plot.allele.frequency <- function(genome,freqs,locus.position){
  allele.frequency <- get.allele.frequency(genome,freqs,locus.position)
  max.freq <- max(allele.frequency)
  allele.number <- get.nb.alleles.per.locus(genome)[locus.position]
  plot(allele.frequency[1,],type="l",ylim=c(0,1.2*max.freq),col=1,xlab = "Generation",ylab="frequency")
  for(allele in 2:allele.number){
    lines(allele.frequency[allele,],col=allele)
  }
  legend("topright",legend=as.character(1:allele.number),lty = rep(1,allele.number),col=1:allele.number)
}

get.allele.frequency <- function(genome,freqs,locus.position)
{
  nb.generation <- ncol(freqs)
  allele.number <- get.nb.alleles.per.locus(genome)[locus.position]
  print(allele.number)
  allele.frequency <- matrix(0,ncol = nb.generation,nrow = allele.number)
  for (generation in 1:nb.generation){
    allele.frequency[,generation] <- get.allele.frequency.single.generation(genome,freqs[,generation],locus.position)
  }
  return(allele.frequency)
}

get.allele.frequency.single.generation <- function(genome,freqs,locus.position){
  allele.number <- get.nb.alleles.per.locus(genome)[locus.position]
  allele.frequency <- sapply(1:allele.number,get.single.allele.frequency.single.generation,genome = genome, freqs = freqs, locus.position = locus.position)
  return(allele.frequency)
}

get.single.allele.frequency.single.generation <- function(allele,genome,freqs,locus.position){
  all.haplotype <- genome@all.haplotype
  haplotype.frequency <- get.haplotype.frequency.single.generation(genome,freqs)
  matching.haplotype <- which(all.haplotype[,locus.position] == allele)
  allele.frequency <- sum(haplotype.frequency[matching.haplotype])
  return(allele.frequency)
}

get.frequency.from.one.allele.frequency <- function(genome,locus,allele,allele.frequency){
  nb.genotype <- get.nb.genotype(genome)
  all.matching.genotype <- get.genotype.with.given.allele(genome,locus,allele)
  matching.homozygothe <- all.matching.genotype[duplicated(all.matching.genotype)]
  all.matching.genotype.not.double <- unique(all.matching.genotype)
  matching.heterozygothe <- setdiff(all.matching.genotype.not.double,matching.homozygothe)
  no.matching <- setdiff(1:nb.genotype,union(matching.homozygothe,matching.heterozygothe))
  nb.matching.genotype <- length(all.matching.genotype)
  frequency <- rep(0, nb.genotype)
  frequency[matching.homozygothe] <- 2*allele.frequency/length(all.matching.genotype)
  frequency[matching.heterozygothe] <- allele.frequency/length(all.matching.genotype)+(1-allele.frequency)/(2*nb.genotype - length(all.matching.genotype))
  frequency[no.matching] <- 2*(1-allele.frequency)/(2*nb.genotype - length(all.matching.genotype))
  return(frequency)
}

get.genotype.with.given.allele <- function(genome,locus,allele){
   all.matching.haplotype <- get.haplotype.with.given.allele(genome,locus,allele)
   all.matching.genotype <- c()
   all.genotype <- genome@all.genotype
   for(haplotype in all.matching.haplotype){
     position.first.haplotype.match <- which(haplotype==all.genotype[,1])
     position.second.haplotype.match <- which(haplotype==all.genotype[,2])
     all.matching.genotype <- c(position.first.haplotype.match,position.second.haplotype.match,all.matching.genotype)
   }
   return(all.matching.genotype)
}

get.haplotype.with.given.allele <- function(genome,locus,allele){
  all.haplotype <- genome@all.haplotype
  all.allele <- all.haplotype[,locus]
  matching.haplotype <- which(all.allele==allele)
  return(matching.haplotype)
}
