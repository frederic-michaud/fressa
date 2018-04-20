
#' get the recombination modifier of a single genotype
#'
get.recombination.modifier.from.genotype <- function(genotype, genome){
  all.genotype <- genome@all.genotype
  all.haplotype <- genome@all.haplotype
  recombination <- 1
  haplotype1 <- get.haplotype.from.index(all.genotype[genotype,1],get.nb.alleles.per.locus(genome))
  haplotype2 <- get.haplotype.from.index(all.genotype[genotype,2],get.nb.alleles.per.locus(genome))
  for (locus in 1:get.nb.locus(genome)){
    position <- where.is.locus(c(haplotype1[locus],haplotype2[locus]),build.genotype.from.locus(genome,locus))
    locus.recombination <- genome@locus[[locus]]$recombination.modifier[position]
    recombination <- recombination*locus.recombination
  }
  return(recombination)
}

#' build the recombination factor of all genotype
#'

build.all.recombination.modifier <- function(genome){
  recombination.factor <- sapply(1:get.nb.genotype(genome), get.recombination.modifier.from.genotype,genome = genome)
  return(recombination.factor)
}

#' get the recombination factor of all genotype
#'
get.all.recombination.modifier <- function(genome){
  return(genome@all.recombination.modifier)
}

