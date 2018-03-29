#' give the position of the sd locus
#'
#' This function return the position of the sexual determining
#' locus, i.e. the locus containing the sex determination, in the genome.
#' If it is not found, return an error.
#' If more than one sd is found, return an error

get.id.sd.locus <- function(genome){
  sd.locus <- 0
  nb.locus <- get.nb.locus(genome)
  for(n.locus in 1:nb.locus){
    is.in <- ("sd" %in% colnames(genome[[n.locus]]))
    if(sd.locus > 0 & is.in) stop("too many sex loci")
    if(is.in) sd.locus <- n.locus

  }
  if(sd.locus ==0) stop("no sexual locus found to determine the sex")
  return(sd.locus)
}

#' Give the list of maleness in a population
get.all.maleness <- function(genome){
  all.haplotype <- build.all.haplotype(genome)
  all.genotype <- build.all.genotype(genome)
  male = c()
  for(genotype.index in 1:get.nb.genotype(genome)){
    male <- c(male,get.maleness(genotype.index,genome))
  }
  return(male)
}

#' Determine the maleness of a genotype
#'
#' Gives the proportion of individuals whith the given genotype
#' which are male
get.maleness <- function(genotype.index, genome){
  all.haplotype <- build.all.haplotype(genome)
  all.genotype <- build.all.genotype(genome)
  sd.locus <- get.id.sd.locus(genome)
  haplotype1 = all.haplotype[all.genotype[genotype.index,1],]
  haplotype2 = all.haplotype[all.genotype[genotype.index,2],]
  my_locus.sd <- genome[[sd.locus]]$sd
  position <- where.is.locus(c(haplotype1[sd.locus],haplotype2[sd.locus]),build.genotype.from.locus(genome,sd.locus))
  return(my_locus.sd[position])
}


#' give the list of male in a population
#'
#' This function return the index of all genotype that are,
#' at least partially, male
#' the number of sd represent of much male is an individual.
#' if it's one, the individual is fully male, if 0 fully female
get.male <- function(genome){
  all.haplotype <- build.all.haplotype(genome)
  all.genotype <- build.all.genotype(genome)
  male = c()
  for(genotype in 1:get.nb.genotype(genome)){
    if(is.male(genotype,genome)) male <- c(male,genotype)
  }
  return(male)
}

#' Determine if a genotype is male
#'
#' This function return true if the genotype with index
#' genotype.index might generate a male and false otherwise.
is.male <- function(genotype.index, genome){
  all.haplotype <- build.all.haplotype(genome)
  all.genotype <- build.all.genotype(genome)
  sd.locus <- get.id.sd.locus(genome)
  haplotype1 = all.haplotype[all.genotype[genotype.index,1],]
  haplotype2 = all.haplotype[all.genotype[genotype.index,2],]
  my_locus.sd <- genome[[sd.locus]]$sd
  position <- where.is.locus(c(haplotype1[sd.locus],haplotype2[sd.locus]),build.genotype.from.locus(genome,sd.locus))
  return(my_locus.sd[position] > 0)
}
#' give the list of female in a population
#'
#' This function return the index of all genotype that are,
#' at least partially, female
get.female <- function(genome){
  all.haplotype <- build.all.haplotype(genome)
  all.genotype <- build.all.genotype(genome)
  female = c()
  for(genotype in 1:get.nb.genotype(genome)){
    if(is.female(genotype, genome)) female <- c(female,genotype)
  }
  return(female)
}

#' Determine if a genotype is female
#'
#' This function return true if the genotype with index
#' genotype.index might generate a female, and false otherwise.
is.female <- function(genotype.index, genome){
  all.haplotype <- build.all.haplotype(genome)
  all.genotype <- build.all.genotype(genome)
  sd.locus <- get.id.sd.locus(genome)
  haplotype1 = all.haplotype[all.genotype[genotype.index,1],]
  haplotype2 = all.haplotype[all.genotype[genotype.index,2],]
  my_locus.sd <- genome[[sd.locus]]$sd
  position <- where.is.locus(c(haplotype1[sd.locus],haplotype2[sd.locus]),build.genotype.from.locus(genome,sd.locus))
  return(my_locus.sd[position] < 1)
}


#' Give the list of femaleness in a population
get.all.femaleness <- function(genome){
  return(1-get.all.maleness(genome))
}

#' Determine the femaleness of a genotype
#'
#' Gives the proportion of individuals whith the given genotype
#' which are female
get.femaleness <- function(genotype.index, genome){
  return(1-get.maleness(genotype.index, genome))
}
