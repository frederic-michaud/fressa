
#' get fitness of a male
#'
#' get the fitness of a male from the index value of its genotype
get.fitness.from.genotype.male <- function(genotype, genome){
  all.genotype <- genome@all.genotype
  all.haplotype <- genome@all.haplotype
  fitness <- 1
  haplotype1 <- get.haplotype.from.index(all.genotype[genotype,1],get.nb.alleles.per.locus(genome))
  haplotype2 <- get.haplotype.from.index(all.genotype[genotype,2],get.nb.alleles.per.locus(genome))
  for (locus in 1:get.nb.locus(genome)){
    position <- where.is.locus(c(haplotype1[locus],haplotype2[locus]),build.genotype.from.locus(genome,locus))
    locus.fitness <- genome@locus[[locus]]$fitness.male[position]
    fitness <- fitness*locus.fitness
  }
  return(fitness)
}



#' get fitness of a female
#'
#' get the fitness of a female from the index value of its genotype
get.fitness.from.genotype.female <- function(genotype, genome){
  all.genotype <- genome@all.genotype
  all.haplotype <- genome@all.haplotype
  fitness <- 1
  haplotype1 <- get.haplotype.from.index(all.genotype[genotype,1],get.nb.alleles.per.locus(genome))
  haplotype2 <- get.haplotype.from.index(all.genotype[genotype,2],get.nb.alleles.per.locus(genome))
  for (locus in 1:get.nb.locus(genome)){
    position <- where.is.locus(c(haplotype1[locus],haplotype2[locus]),build.genotype.from.locus(genome,locus))
    locus.fitness <- genome@locus[[locus]]$fitness.female[position]
    fitness <- fitness*locus.fitness
  }
  return(fitness)
}

build.all.fitness.female <- function(genome){
  female.genotype <- sapply(1:get.nb.genotype(genome), get.fitness.from.genotype.female,genome = genome)
  return(female.genotype)
}

build.all.fitness.male <- function(genome){
  male.genotype <- sapply(1:get.nb.genotype(genome), get.fitness.from.genotype.male,genome = genome)
  return(male.genotype)
}

get.all.fitness.female <- function(genome){
  return(genome@all.fitness.female)
}

get.all.fitness.male <- function(genome){
  return(genome@all.fitness.male)
}
