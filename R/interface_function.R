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
  palette = get.palette(get.nb.haplotype(genome))
  plot(haplotype.frequency[1,],type="l",ylim=c(0,1.2*max.freq),col=palette[1],xlab = "Generation",ylab="frequency")
  for(haplotype in 2:get.nb.haplotype(genome)){
    lines(haplotype.frequency[haplotype,],col=palette[haplotype])
  }
  legend("topright",legend=get.haplotype.names(genome),lty = rep(1,get.nb.haplotype(genome)),col=palette)
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
  palette <- get.palette(allele.number)
  plot(allele.frequency[1,],type="l",ylim=c(0,1.2*max.freq),col=palette[1],xlab = "Generation",ylab="frequency")
  for(allele in 2:allele.number){
    lines(allele.frequency[allele,],col=palette[allele])
  }
  legend("topright",legend=get.allele.name(genome,locus.position),lty = rep(1,allele.number),col=palette)
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
  matching.haplotype <- get.haplotype.with.given.allele(genome,locus,allele)

  male.haplotype <- get.haplotype.male(genome)
  male.matching.haplotype = intersect(matching.haplotype,male.haplotype)
  female.haplotype <- get.haplotype.female(genome)
  female.matching.haplotype = intersect(matching.haplotype,female.haplotype)

  nb.haplotype <- get.nb.haplotype(genome)
  nb.male.haplotype <- length(male.haplotype)
  nb.female.haplotype <- length(female.haplotype)
  nb.male.matching.haplotype <- length(male.matching.haplotype)
  nb.female.matching.haplotype <- length(female.matching.haplotype)


  haplotype.matching.male.frequency <- 0.5*allele.frequency/nb.male.matching.haplotype
  haplotype.matching.female.frequency <- 0.5*allele.frequency/nb.female.matching.haplotype

  if(nb.male.matching.haplotype > 0){
    haplotype.no.machting.male.frequency <- (1-0.5*allele.frequency)/(nb.male.haplotype - nb.male.matching.haplotype)

  }
  else{
    haplotype.no.machting.male.frequency <- 1/nb.male.haplotype
  }

  if(nb.female.matching.haplotype > 0){
    haplotype.no.matching.female.frequency <- (1-0.5*allele.frequency)/(nb.female.haplotype - nb.female.matching.haplotype)
  }
  else{
    haplotype.no.matching.female.frequency <- 1/nb.female.haplotype
  }

  male.haplotype.frequency <- rep(0,nb.haplotype)
  female.haplotype.frequency <- rep(0,nb.haplotype)

  male.haplotype.frequency[male.haplotype] <- haplotype.no.machting.male.frequency
  female.haplotype.frequency[female.haplotype] <- haplotype.no.matching.female.frequency


  male.haplotype.frequency[male.matching.haplotype] <- haplotype.matching.male.frequency
  female.haplotype.frequency[female.matching.haplotype] <- haplotype.matching.female.frequency

  genotype.frequency.as.matrix <- outer(male.haplotype.frequency,female.haplotype.frequency)
  genotype.frequency.as.matrix <-  genotype.frequency.as.matrix +t(genotype.frequency.as.matrix)
  diag(genotype.frequency.as.matrix) <- diag(genotype.frequency.as.matrix)/2
  frequency <- genotype.frequency.as.matrix[genome@all.genotype]
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


