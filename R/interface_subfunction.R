#' return names for the haplotype if a name is present on a locus
get.haplotype.names.from.allele.names <- function(genome){
  all.haplotype <- genome@all.haplotype
  haplotype.names <- c()
  for(haplotype in 1:get.nb.haplotype(genome)){
    haplotype.name <- c()
    haplotype <- all.haplotype[haplotype,]
    for(locus in 1:get.nb.locus(genome)){
      haplotype.name <- paste(haplotype.name, genome@locus[[locus]]@allele.name[haplotype[locus]],sep = "",collapse = "")
    }
    haplotype.names <- c(haplotype.names,haplotype.name)
  }
  return(haplotype.names)
}

#' return names for the haplotype if no name is provided
get.haplotype.names.from.allele.number <- function(genome){
  all.haplotype <- genome@all.haplotype
  haplotype.names <- apply(all.haplotype,1,paste,collapse="")
  return(haplotype.names)
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

#get the marignal fitness of all the haplotype for a given generation

get.marginal.haplotype.fitness.single.generation <- function(genome,frequency){
  nb.haplotype <- get.nb.haplotype(genome)
  marginal.fitnesses <- sapply(1:nb.haplotype,function(iter) get.marginal.haplotype.fitness.single.generation.single.haplotype(genome,frequency, iter))
  return(marginal.fitnesses)
}

#get the marignal fitness of one of the haplotype for a given generation

get.marginal.haplotype.fitness.single.generation.single.haplotype <- function(genome,frequency,haplotype){
  matching.genotype <- get.genotype.with.given.haplotype(genome,haplotype)

  all.fitness.in.male <- get.all.fitness.male(genome)
  maleness <- get.all.maleness(genome)
  marginal.fitness.male <- sum(frequency[matching.genotype]*maleness[matching.genotype]*all.fitness.in.male[matching.genotype])

  all.fitness.in.female <- get.all.fitness.female(genome)
  femaleness <- get.all.femaleness(genome)
  marginal.fitness.female<- sum(frequency[matching.genotype]*femaleness[matching.genotype]*all.fitness.in.female[matching.genotype])

  overall.marginal.fitness <- (marginal.fitness.male + marginal.fitness.female)/sum(frequency[matching.genotype])

  return(overall.marginal.fitness)
}

#get the frequency of an haplotype in one generation

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

#get the frequency of all allele in one generation for a given locus


get.allele.frequency.single.generation <- function(genome,freqs,locus.position){
  allele.number <- get.nb.alleles.per.locus(genome)[locus.position]
  allele.frequency <- sapply(1:allele.number,get.single.allele.frequency.single.generation,genome = genome, freqs = freqs, locus.position = locus.position)
  return(allele.frequency)
}

#get the frequency of one allele in one generation

get.single.allele.frequency.single.generation <- function(allele,genome,freqs,locus.position){
  all.haplotype <- genome@all.haplotype
  haplotype.frequency <- get.haplotype.frequency.single.generation(genome,freqs)
  matching.haplotype <- which(all.haplotype[,locus.position] == allele)
  allele.frequency <- sum(haplotype.frequency[matching.haplotype])
  return(allele.frequency)
}
