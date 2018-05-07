#' return names for the gamete if a name is present on a locus
get.gamete.names.from.allele.names <- function(genome){
  all.gamete <- genome@all.gamete
  gamete.names <- c()
  for(gamete in 1:get.nb.gamete(genome)){
    gamete.name <- c()
    gamete <- all.gamete[gamete,]
    for(locus in 1:get.nb.locus(genome)){
      gamete.name <- paste(gamete.name, genome@locus[[locus]]@allele.name[gamete[locus]],sep = "",collapse = "")
    }
    gamete.names <- c(gamete.names,gamete.name)
  }
  return(gamete.names)
}

#' return names for the gamete if no name is provided
get.gamete.names.from.allele.number <- function(genome){
  all.gamete <- genome@all.gamete
  gamete.names <- apply(all.gamete,1,paste,collapse="")
  return(gamete.names)
}



#get the names of the allele for a given locus
get.allele.name <- function(genome,locus){
  if (length(genome@locus[[locus]]@allele.name) > 0) names <- genome@locus[[locus]]@allele.name
  else names <-  as.character(1:get.nb.alleles.per.locus(genome)[locus])
}

#get the marginal fitness of all the allele of a locus for a single generation
get.marginal.allele.fitness.single.generation <- function(genome,frequency,locus){
  nb.allele <- get.nb.alleles.per.locus(genome)[locus]
  marginal.fitnesses <- sapply(1:nb.allele,function(iter) get.marginal.allele.fitness.single.generation.single.allele(genome,frequency, locus,iter))
  return(marginal.fitnesses)
}

#get the marginal fitness of one allele for a single generation

get.marginal.allele.fitness.single.generation.single.allele <- function(genome,frequency,locus,allele){
  matching.genotype <- get.genotype.with.given.allele(genome,locus,allele)

  all.fitness.in.male <- get.all.fitness.male(genome)
  maleness <- get.all.maleness(genome)
  marginal.fitness.male <- sum(frequency[matching.genotype]*maleness[matching.genotype]*all.fitness.in.male[matching.genotype])

  all.fitness.in.female <- get.all.fitness.female(genome)
  femaleness <- get.all.femaleness(genome)
  marginal.fitness.female<- sum(frequency[matching.genotype]*femaleness[matching.genotype]*all.fitness.in.female[matching.genotype])

  overall.marginal.fitness <- (marginal.fitness.male + marginal.fitness.female)/sum(frequency[matching.genotype])

  return(overall.marginal.fitness)
}

#get the marignal fitness of all the gamete for a given generation

get.marginal.gamete.fitness.single.generation <- function(genome,frequency){
  nb.gamete <- get.nb.gamete(genome)
  marginal.fitnesses <- sapply(1:nb.gamete,function(iter) get.marginal.gamete.fitness.single.generation.single.gamete(genome,frequency, iter))
  return(marginal.fitnesses)
}

#get the marignal fitness of one of the gamete for a given generation

get.marginal.gamete.fitness.single.generation.single.gamete <- function(genome,frequency,gamete){
  matching.genotype <- get.genotype.with.given.gamete(genome,gamete)

  all.fitness.in.male <- get.all.fitness.male(genome)
  maleness <- get.all.maleness(genome)
  marginal.fitness.male <- sum(frequency[matching.genotype]*maleness[matching.genotype]*all.fitness.in.male[matching.genotype])

  all.fitness.in.female <- get.all.fitness.female(genome)
  femaleness <- get.all.femaleness(genome)
  marginal.fitness.female<- sum(frequency[matching.genotype]*femaleness[matching.genotype]*all.fitness.in.female[matching.genotype])

  overall.marginal.fitness <- (marginal.fitness.male + marginal.fitness.female)/sum(frequency[matching.genotype])

  return(overall.marginal.fitness)
}

#get the frequency of an gamete in one generation

get.gamete.frequency.single.generation <- function(genome,freqs)
{
  nb.gamete <- get.nb.gamete(genome)
  nb.genotype <- get.nb.genotype(genome)
  all.genotype <- genome@all.genotype
  #frequency is a matrix with the frequency of the various genome
  #row index indicate first gamete and column index indicate second gamete
  frequency <- matrix(0,nrow = nb.gamete,ncol=nb.gamete)
  for (genotype in 1:nb.genotype){
    frequency[all.genotype[genotype,1],all.genotype[genotype,2]] <- freqs[genotype]
  }
  sum.column <- colSums(frequency)/2
  sum.row <- rowSums(frequency)/2
  return(sum.column + sum.row)
}

#get the frequency of all allele in one generation for a given locus


get.allele.frequency.single.generation <- function(genome,freqs,locus.position){
  allele.number <- get.nb.alleles.per.locus(genome)[locus.position]
  allele.frequency <- sapply(1:allele.number,get.single.allele.frequency.single.generation,genome = genome, freqs = freqs, locus.position = locus.position)
  return(allele.frequency)
}

#get the frequency of one allele in one generation

get.single.allele.frequency.single.generation <- function(allele,genome,freqs,locus.position){
  all.gamete <- genome@all.gamete
  gamete.frequency <- get.gamete.frequency.single.generation(genome,freqs)
  matching.gamete <- which(all.gamete[,locus.position] == allele)
  allele.frequency <- sum(gamete.frequency[matching.gamete])
  return(allele.frequency)
}

#has the simulation converged?

is.converged <- function(freqs,generation,min.generations,max.generations,criteria){
  #If we are still in the warmup, we wait
  if(generation < min.generations) return(FALSE)

  last.slope <- sum(abs(freqs[,generation] - freqs[,generation-1]))
  pre.last.slope <- sum(abs(freqs[,generation-1] - freqs[,generation-2]))
  #if the criteria is increasing, we should better wait
  if(last.slope > pre.last.slope) return(FALSE)

  if(generation > max.generations & max.generations > 0){
    warning("The maximum number of generation has been reached without reaching convergence")
    return(TRUE)
  }

  is.criteria.met <- (last.slope < criteria)
  return(is.criteria.met)
}
