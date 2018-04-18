#' Get the genotype of all child from their parents
#'
#' This function returns a 2x4 matrix containing all the possible
#' genotype (on each line) of the child coming from two differents
#' parents.
get.childs.genotype.from.parent.genotype <- function(male.genotype,female.genotype,genome){
  haplotype.male.1 <- male.genotype[1]
  haplotype.male.2 <- male.genotype[2]
  haplotype.female.1 <- female.genotype[1]
  haplotype.female.2 <- female.genotype[2]
  a1 = c(haplotype.male.1,haplotype.female.1)
  a2 = c(haplotype.male.1,haplotype.female.2)
  a3 = c(haplotype.male.2,haplotype.female.1)
  a4 = c(haplotype.male.2,haplotype.female.2)
  return(matrix(c(a1,a2,a3,a4),ncol=2,byrow = T))
}


#' Get the index of all child from their parents index
#'
#' This function returns the index of all child, i.e. the position
#' of this individual in the all.genotype matrix.
get.childs.index.from.parent.index <- function(male.genotype,female.genotype,genome){
  all.genotype <- genome@all.genotype
  childs <- get.childs.genotype.from.parent.genotype(all.genotype[male.genotype,],all.genotype[female.genotype,])
  a1 <- get.genotype.index.from.haplotypes.index(childs[1,],all.genotype)
  a2 <- get.genotype.index.from.haplotypes.index(childs[2,],all.genotype)
  a3 <- get.genotype.index.from.haplotypes.index(childs[3,],all.genotype)
  a4 <- get.genotype.index.from.haplotypes.index(childs[4,],all.genotype)
  return(c(a1,a2,a3,a4))
}

#' Get the frequency of new born as a function of adults

simulate.frequency <- function(genome,initial.frequency){
  female.gamete.frequency <- get.female.gamete.frequency(genome,initial.frequency)
  male.gamete.frequency <- get.male.gamete.frequency(genome,initial.frequency)
  genotype.frequency.as.matrix <- outer(male.gamete.frequency,female.gamete.frequency)
  genotype.frequency.as.matrix <-  genotype.frequency.as.matrix +t(genotype.frequency.as.matrix)
  diag(genotype.frequency.as.matrix) <- diag(genotype.frequency.as.matrix)/2
  frequency <- genotype.frequency.as.matrix[genome@all.genotype]
  return(frequency)
}



