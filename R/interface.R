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
#' genome = create.genome(locus1,locus2)
#' freqs <- compute.frequency.evolution(genome)
#' @export


compute.frequency.evolution <- function(genome,initial.frequency = NULL,generations = 25)
{
  nb.genotype <- get.nb.genotype(genome)
  if(is.null(initial.frequency)){
    initial.frequency <- rep(1/nb.genotype,nb.genotype)
  }
  if(length(initial.frequency)!=nb.genotype) stop("The given initial frequency is not the correct size")
  freqs <- matrix(0,nrow = nb.genotype,ncol = generations)
  freqs[,1] <- initial.frequency
  for(generation in 2:generations)
  {
    freqs[,generation] <- simulate.frequency(genome,freqs[ ,generation-1])
  }
  return(freqs)
}

#' Return the haplotype frequency in a population
#'
#' given a matrix of frequency returned by the function `compute.frequency.evolution`
#' and the associated genome, return a matrix containing the frequency of each haplotype
#' through time. Each row contains a genotype while each column contains a generation.
#'
#' @param genome A S4 object of type genome
#' @param freqs a matrix of frequency as returned by the function `compute.frequency.evolution`
#' @examples
#' locus1 = create.locus(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1))
#' locus2 = create.locus(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#' freqs <- compute.frequency.evolution(genome)
#' freqs.haplotype <- get.haplotype.frequency(genome, freqs)
#' @export

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

#' Return the allele frequency for a given locus
#'
#' given a matrix of frequency returned by the function `compute.frequency.evolution`
#' and the associated genome, return a matrix containing the frequency of each allele
#' of a given locuc through time. Each row contains an allele while each column
#' contains a generation.
#'
#' @param genome A S4 object of type genome
#' @param freqs a matrix of frequency as returned by the function `compute.frequency.evolution`
#' @examples
#' locus1 = create.locus(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1))
#' locus2 = create.locus(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#' freqs <- compute.frequency.evolution(genome)
#' freqs.allele <- get.allele.frequency(genome, freqs)
#' @export

get.allele.frequency <- function(genome,freqs,locus.position)
{
  nb.generation <- ncol(freqs)
  allele.number <- get.nb.alleles.per.locus(genome)[locus.position]
  allele.frequency <- matrix(0,ncol = nb.generation,nrow = allele.number)
  for (generation in 1:nb.generation){
    allele.frequency[,generation] <- get.allele.frequency.single.generation(genome,freqs[,generation],locus.position)
  }
  return(allele.frequency)
}

#' This functions allows to generate genotypic frequency
#' by specifying the initial frequency of one allele.
#'
#'This function is usefull mainly to start a simulation with one allele being
#'rare. If no input frequency is given, it will be assumed that all genotype have
#'equal frequency. However, it might be usefull to start with one allele being
#'very rare. This function allows to do so by automatically computing a
#'initial set of frequency. To do so, it produce a pool of gamete where the given allele
#'appears with the specified allele, and then perform random mating (without taking into account fitness).
#'This ensure that the allele is present at the correct level, and that there is no bias in the
#'sex-ratio.
#'
#' @param genome A S4 object of type genome
#' @param locus the locus at which the allele should have a given frequency
#' @param allele the allele which has a given frequency
#' @param allele.frequency the initial frequency of the allele.
#' @examples
#' locus1 = create.locus(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1))
#' locus2 = create.locus(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#' initial.frequency <- get.frequency.from.one.allele.frequency(genome,2,1,0.01)
#' freqs <- compute.frequency.evolution(genome,initial.frequency)
#' freqs.allele <- get.allele.frequency(genome, freqs)
#' @export

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

#' get the marginal fitness of all haplotype in the population
#'
#' given a matrix of frequency returned by the function `compute.frequency.evolution`
#' and the associated genome, return a matrix containing the evolution of the marginal
#' fitness. The marginal fitness is defined as the mean fitness of
#' individual carrying this haplotype weighted by the frequency of those individuals.
#'
#' @param genome A S4 object of type genome
#' @param freqs a matrix of frequency as returned by the function `compute.frequency.evolution`
#' @examples
#' locus1 = create.locus(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1))
#' locus2 = create.locus(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#' freqs <- compute.frequency.evolution(genome)
#' get.haplotype.marginal.fitness(genome, freqs)
#' @export


get.marginal.haplotype.fitness <- function(genome,freqs)
{
  nb.generation <- ncol(freqs)
  haplotype.number <- get.nb.haplotype(genome)
  haplotype.marginal.fitness <- matrix(0,ncol = nb.generation,nrow = haplotype.number)
  for (generation in 1:nb.generation){
    haplotype.marginal.fitness[,generation] <- get.marginal.haplotype.fitness.single.generation(genome,freqs[,generation])
  }
  return(haplotype.marginal.fitness)
}

#' get the marginal fitness of all allele from one locus
#'
#' given a matrix of frequency returned by the function `compute.frequency.evolution`
#' and the associated genome, return a matrix containing the evolution of the marginal
#' fitness of all allele. The marginal fitness is defined as the mean fitness of
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
#' get.haplotype.marginal.fitness(genome, freqs)
#' @export

get.marginal.allele.fitness <- function(genome,freqs,locus)
{
  nb.generation <- ncol(freqs)
  allele.number <- get.nb.alleles.per.locus(genome)[locus]
  allele.marginal.fitness <- matrix(0,ncol = nb.generation,nrow = allele.number)
  for (generation in 1:nb.generation){
    allele.marginal.fitness[,generation] <- get.marginal.allele.fitness.single.generation(genome,freqs[,generation],locus)
  }
  return(allele.marginal.fitness)
}

#' get the name of all genotype present in the population
#'
#' This function returns a name for all genotype present in the population
#' This is useful for plotting result but also to know in which order the genotype
#' are store in genome.  This is useful for example to specify the initial frequency
#'
#' Notice that if allele.name is specified, the name will contain this name and is
#' therefore much easier to read that if it's not present, where the name of the allele
#' is just their number.
#'
#' @param genome A S4 object of type genome

#' @examples
#' locus1 = create.locus(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1),allele.name = c("x","y"))
#' locus2 = create.locus(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1),allele.name = c("F","M"))
#' genome = create.genome(locus=list(locus1,locus2))
#' get.genotype.names(genome)
#' @export

get.genotype.names <- function(genome){
  all.genotype <- genome@all.genotype
  haplotype.names <- get.haplotype.names(genome)
  genotype.names=c()
  for(genotype in 1:get.nb.genotype(genome)){
    haplotype1.name <- haplotype.names[all.genotype[genotype,1]]
    haplotype2.name <- haplotype.names[all.genotype[genotype,2]]
    genotype.name <-paste(haplotype1.name,haplotype2.name, sep="|")
    genotype.names <- c(genotype.names,genotype.name)
  }
  return(genotype.names)
}

#' get the name of all haplotype present in the population
#'
#' This function returns a name for all haplotype present in the population
#' This is useful mainly for plotting result but also to know in which order the haplotype
#' are store in genome.
#'
#' Notice that if allele.name is specified, the name will contain this name and is
#' therefore much easier to read that if it's not present, where the name of the allele
#' is just their number.
#'
#' @param genome A S4 object of type genome

#' @examples
#' locus1 = create.locus(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1),allele.name = c("x","y"))
#' locus2 = create.locus(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1),allele.name = c("F","M"))
#' genome = create.genome(locus=list(locus1,locus2))
#' get.haplotype.names(genome)
#' @export

get.haplotype.names <- function(genome){
  if(length(genome@locus[[1]]@allele.name) == 0) haplotype.names <- get.haplotype.names.from.allele.number(genome)
  else haplotype.names <- get.haplotype.names.from.allele.names(genome)
  return(haplotype.names)
}

