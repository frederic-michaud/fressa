#' return names for the haplotype
get.haplotype.names <- function(genome){
  if(length(genome@locus[[1]]@allele.name) == 0) haplotype.names <- get.haplotype.names.without.names(genome)
  else haplotype.names <- get.haplotype.names.with.names(genome)
  return(haplotype.names)
}

#' return names for the haplotype if a name is present on a locus
get.haplotype.names.with.names <- function(genome){
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
get.haplotype.names.without.names <- function(genome){
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

#get the names of the allele
get.allele.name <- function(genome,locus){
  if (length(genome@locus[[locus]]@allele.name) > 0) names <- genome@locus[[locus]]@allele.name
  else names <-  as.character(1:get.nb.alleles.per.locus(genome)[locus])
}

#get a palette of ncolor colors as much different as possible from each other
get.palette <- function(n.color){
  palette <- rainbow(n.color)[sample(1:n.color)]
  return(palette)
}
