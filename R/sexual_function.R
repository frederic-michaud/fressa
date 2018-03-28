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

get.male <- function(genome,all.haplotype,all.genotype){
  sd.locus <- get.id.sd.locus(genome)
  male = c()
  for(genotype in 1:get.nb.genotype(genome)){
    if(is.male(genotype,sd.locus,all.haplotype,all.genotype,genome)) male <- c(male,genotype)
  }
  return(male)
}

is.male <- function(genotype,sd.locus,all.haplotype,all.genotype,genome){
  haplotype1 = all.haplotype[all.genotype[genotype,1]]
  haplotype2 = all.haplotype[all.genotype[genotype,2]]
  my_locus.sd <- genome[[sd.locus]]$sd
  position <- where.is.locus(c(haplotype1[sd.locus],haplotype2[sd.locus]),build.genotype.from.locus(genome,sd.locus))
  return(my_locus.sd[position] == 1)
}

get.female <- function(genome,all.haplotype,all.genotype){
  sd.locus <- get.id.sd.locus(genome)
  female = c()
  for(genotype in 1:get.nb.genotype(genome)){
    if(is.female(genotype,sd.locus,all.haplotype,all.genotype,genome)) female <- c(female,genotype)
  }
  return(female)
}

is.female <- function(genotype,sd.locus,all.haplotype,all.genotype,genome){
  haplotype1 = all.haplotype[all.genotype[genotype,1]]
  haplotype2 = all.haplotype[all.genotype[genotype,2]]
  my_locus.sd <- genome[[sd.locus]]$sd
  position <- where.is.locus(c(haplotype1[sd.locus],haplotype2[sd.locus]),build.genotype.from.locus(genome,sd.locus))
  return(my_locus.sd[position] == 0)
}


