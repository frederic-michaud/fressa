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
#' locus1 = data.frame(allele1=c(1,1),allele2 = c(1,2),sd = c(0,1),fitness.male=c(1,1),fitness.female=c(1,1))
#' locus2 = data.frame(allele1=  c(1,1,2),allele2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
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


#' Compute the evolution of the frequency of the genotypes until convergence is reach
#'
#' Given a genome, this function simulate the evolution of the frequency
#' of all the possible genotype that exist until convergence has been reach.
#'
#' Convergence is evaluated in the following way. First, a few generation are perform (warmup) without measuring the criteria
#' since they might be oscillation due to male/female proportion effect. Then, at each generation, the total
#' slope $\abs(\vec({_t}-\vec{x_{t-1}})$ and the total curvature $\abs(\vec({_t}-\vec{x_{t-1}})$
#' @return A matrix containing the frequencies of each genotype at each generation
#' @param genome A S4 object of the type genome
#' @param initial.frequency The initial frequency of the various genotype. If NULL is given
#' the initial frequencies will all be set to the same value
#' @param min.generations The minimum number of generation to be computed
#' @param criteria If the sum of the slope of the evolution of the genotype frequency is under this number, the simulation stop (if the other criteria are also met)
#' @examples
#' locus1 = data.frame(allele1=c(1,1),allele2 = c(1,2),sd = c(0,1),fitness.male=c(1,1),fitness.female=c(1,1))
#' locus2 = data.frame(allele1=  c(1,1,2),allele2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus1,locus2)
#' freqs <- compute.frequency.evolution.until.convergence(genome)
#' @export


compute.frequency.evolution.until.convergence <- function(genome,
                                                          initial.frequency = NULL,
                                                          min.generations = 25,
                                                          criteria = 1.e-8
                                                          )
{
  nb.genotype <- get.nb.genotype(genome)
  if(is.null(initial.frequency)){
    initial.frequency <- rep(1/nb.genotype,nb.genotype)
  }
  if(length(initial.frequency)!=nb.genotype) stop("The given initial frequency is not the correct size")

  freqs <- cbind(initial.frequency,simulate.frequency(genome,initial.frequency)) #we compute the first generation here to be able to use tail
  generation <- 2
  while(!is.converged(freqs,generation,min.generations,criteria))
  {
    generation <- generation +1
    freqs <- cbind(freqs,simulate.frequency(genome,freqs[,generation-1]))

  }
  return(freqs)
}

is.converged <- function(freqs,generation,min.generations,criteria){
  #If we are still in the warmup, we wait
  if(generation < min.generations) return(FALSE)

  last.slope <- sum(abs(freqs[,generation] - freqs[,generation-1]))
  pre.last.slope <- sum(abs(freqs[,generation-1] - freqs[,generation-2]))
  #if the criteria is increasing, we should better wait
  if(last.slope > pre.last.slope) return(FALSE)

  is.criteria.met <- (last.slope < criteria)
  return(is.criteria.met)
}




#' Return the gamete frequency in a population
#'
#' given a matrix of frequency returned by the function `compute.frequency.evolution`
#' and the associated genome, return a matrix containing the frequency of each gamete
#' through time. Each row contains a genotype while each column contains a generation.
#'
#' @param genome A S4 object of type genome
#' @param freqs a matrix of frequency as returned by the function `compute.frequency.evolution`
#' @examples
#' locus1 = create.locus(allele1=c(1,1),allele2 = c(1,2),sd = c(0,1))
#' locus2 = create.locus(allele1=  c(1,1,2),allele2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#' freqs <- compute.frequency.evolution(genome)
#' freqs.gamete <- get.gamete.frequency(genome, freqs)
#' @export

get.gamete.frequency <- function(genome,freqs)
{
  nb.generation <- ncol(freqs)
  gamete.frequency <- matrix(0,ncol = nb.generation,nrow = get.nb.gamete(genome))
  for (generation in 1:nb.generation){
    gamete.frequency[,generation] <- get.gamete.frequency.single.generation(genome,freqs[,generation])
  }
  row.names(gamete.frequency) <- get.gamete.names(genome)
  return(gamete.frequency)
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
#' locus1 = create.locus(allele1=c(1,1),allele2 = c(1,2),sd = c(0,1))
#' locus2 = create.locus(allele1=  c(1,1,2),allele2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
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
#' locus1 = create.locus(allele1=c(1,1),allele2 = c(1,2),sd = c(0,1))
#' locus2 = create.locus(allele1=  c(1,1,2),allele2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#' initial.frequency <- get.frequency.from.one.allele.frequency(genome,2,1,0.01)
#' freqs <- compute.frequency.evolution(genome,initial.frequency)
#' freqs.allele <- get.allele.frequency(genome, freqs)
#' @export

get.frequency.from.one.allele.frequency <- function(genome,locus,allele,allele.frequency){
  matching.gamete <- get.gamete.with.given.allele(genome,locus,allele)

  male.gamete <- get.gamete.male(genome)
  male.matching.gamete = intersect(matching.gamete,male.gamete)
  female.gamete <- get.gamete.female(genome)
  female.matching.gamete = intersect(matching.gamete,female.gamete)

  nb.gamete <- get.nb.gamete(genome)
  nb.male.gamete <- length(male.gamete)
  nb.female.gamete <- length(female.gamete)
  nb.male.matching.gamete <- length(male.matching.gamete)
  nb.female.matching.gamete <- length(female.matching.gamete)


  gamete.matching.male.frequency <- 0.5*allele.frequency/nb.male.matching.gamete
  gamete.matching.female.frequency <- 0.5*allele.frequency/nb.female.matching.gamete

  if(nb.male.matching.gamete > 0){
    gamete.no.machting.male.frequency <- (1-0.5*allele.frequency)/(nb.male.gamete - nb.male.matching.gamete)

  }
  else{
    gamete.no.machting.male.frequency <- 1/nb.male.gamete
  }

  if(nb.female.matching.gamete > 0){
    gamete.no.matching.female.frequency <- (1-0.5*allele.frequency)/(nb.female.gamete - nb.female.matching.gamete)
  }
  else{
    gamete.no.matching.female.frequency <- 1/nb.female.gamete
  }

  male.gamete.frequency <- rep(0,nb.gamete)
  female.gamete.frequency <- rep(0,nb.gamete)

  male.gamete.frequency[male.gamete] <- gamete.no.machting.male.frequency
  female.gamete.frequency[female.gamete] <- gamete.no.matching.female.frequency


  male.gamete.frequency[male.matching.gamete] <- gamete.matching.male.frequency
  female.gamete.frequency[female.matching.gamete] <- gamete.matching.female.frequency

  frequency <- get.frequency.from.gamete.frequency(genome,male.gamete.frequency,female.gamete.frequency)
  return(frequency)
}

#' get the marginal fitness of all gamete in the population
#'
#' given a matrix of frequency returned by the function `compute.frequency.evolution`
#' and the associated genome, return a matrix containing the evolution of the marginal
#' fitness. The marginal fitness is defined as the mean fitness of
#' individual carrying this gamete weighted by the frequency of those individuals.
#'
#' @param genome A S4 object of type genome
#' @param freqs a matrix of frequency as returned by the function `compute.frequency.evolution`
#' @examples
#' locus1 = create.locus(allele1=c(1,1),allele2 = c(1,2),sd = c(0,1))
#' locus2 = create.locus(allele1=  c(1,1,2),allele2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#' freqs <- compute.frequency.evolution(genome)
#' get.gamete.marginal.fitness(genome, freqs)
#' @export


get.marginal.gamete.fitness <- function(genome,freqs)
{
  nb.generation <- ncol(freqs)
  gamete.number <- get.nb.gamete(genome)
  gamete.marginal.fitness <- matrix(0,ncol = nb.generation,nrow = gamete.number)
  for (generation in 1:nb.generation){
    gamete.marginal.fitness[,generation] <- get.marginal.gamete.fitness.single.generation(genome,freqs[,generation])
  }
  return(gamete.marginal.fitness)
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
#' locus1 = create.locus(allele1=c(1,1),allele2 = c(1,2),sd = c(0,1))
#' locus2 = create.locus(allele1=  c(1,1,2),allele2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#' freqs <- compute.frequency.evolution(genome)
#' get.gamete.marginal.fitness(genome, freqs)
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
#' locus1 = create.locus(allele1=c(1,1),allele2 = c(1,2),sd = c(0,1),allele.name = c("x","y"))
#' locus2 = create.locus(allele1=  c(1,1,2),allele2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1),allele.name = c("F","M"))
#' genome = create.genome(locus=list(locus1,locus2))
#' get.genotype.names(genome)
#' @export

get.genotype.names <- function(genome){
  all.genotype <- genome@all.genotype
  gamete.names <- get.gamete.names(genome)
  genotype.names=c()
  for(genotype in 1:get.nb.genotype(genome)){
    gamete1.name <- gamete.names[all.genotype[genotype,1]]
    gamete2.name <- gamete.names[all.genotype[genotype,2]]
    genotype.name <-paste(gamete1.name,gamete2.name, sep="|")
    genotype.names <- c(genotype.names,genotype.name)
  }
  return(genotype.names)
}

#' get the name of all gamete present in the population
#'
#' This function returns a name for all gamete present in the population
#' This is useful mainly for plotting result but also to know in which order the gamete
#' are store in genome.
#'
#' Notice that if allele.name is specified, the name will contain this name and is
#' therefore much easier to read that if it's not present, where the name of the allele
#' is just their number.
#'
#' @param genome A S4 object of type genome

#' @examples
#' locus1 = create.locus(allele1=c(1,1),allele2 = c(1,2),sd = c(0,1),allele.name = c("x","y"))
#' locus2 = create.locus(allele1=  c(1,1,2),allele2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1),allele.name = c("F","M"))
#' genome = create.genome(locus=list(locus1,locus2))
#' get.gamete.names(genome)
#' @export

get.gamete.names <- function(genome){
  if(length(genome@locus[[1]]@allele.name) == 0) gamete.names <- get.gamete.names.from.allele.number(genome)
  else gamete.names <- get.gamete.names.from.allele.names(genome)
  return(gamete.names)
}

#' This functions allows to generate genotypic frequency
#' by specifying the initial frequency og the gamete in male and female
#'
#'This function is usefull mainly to start a simulation with controlled
#'gamete frequency
#'
#' @param genome A S4 object of type genome
#' @param male.gamete.frequency a list of the frequency of the gamete for male (should sum up to one)
#' @param female.gamete.frequency a list of the frequency of the gamete for female (should sum up to one)
#' @export

get.frequency.from.gamete.frequency <- function(genome,male.gamete.frequency,female.gamete.frequency){
  genotype.frequency.as.matrix <- outer(male.gamete.frequency,female.gamete.frequency)
  genotype.frequency.as.matrix <-  genotype.frequency.as.matrix +t(genotype.frequency.as.matrix)
  diag(genotype.frequency.as.matrix) <- diag(genotype.frequency.as.matrix)/2
  frequency <- genotype.frequency.as.matrix[genome@all.genotype]
  return(frequency)
}
