#' return names for the haplotype
get.haplotype.names <- function(genome){
  all.haplotype <- genome@all.haplotype
  haplotype.names <- apply(all.haplotype,1,paste,collapse="")
  return(haplotype.names)
}

#' generate names for the phenotype
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
