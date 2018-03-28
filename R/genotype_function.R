
get.nb.locus <- function(genome){
  return(length(genome))
}

get.nb.genotype <- function(genome)
{
  all.genotype <- build.all.genotype(genome)
  return(dim(all.genotype)[1])
}

#get an haplotype from a number
get.haplotype.from.index <- function(haplotype,nb.alleles){
  representation <- c()
  prod <- 1
  for(n.locus in seq(length(nb.alleles),1,-1)){
    new <- mod(ceiling(haplotype/prod),nb.alleles[n.locus])
    representation <- c(new,representation)
    prod <- prod*nb.alleles[n.locus]
  }
  return(representation)
}

get.nb.alleles.per.locus <- function(genome){
  nb.alleles <- c()
  for (locus in genome){
    alleles <- unique(c(locus$chrom1,locus$chrom2))
    nb.alleles <- c(nb.alleles,length(alleles))
  }
  return(nb.alleles)
}

build.all.haplotype <- function(genome){
  nb.alleles <- get.nb.alleles.per.locus(genome)
  nb.haplotypes <- prod(nb.alleles)
  all.haplotype <- matrix(0,nrow=nb.haplotypes,ncol = length(genome))
  for (haplotype in 1:nb.haplotypes){
    haplotype.representation <- get.haplotype.from.index(haplotype,nb.alleles)
    all.haplotype[haplotype,] <- haplotype.representation
  }
  return(all.haplotype)
}

#For a locus, give all possible genotype, like xx,xy.
build.genotype.from.locus <- function(genome,locus){
   all.config <- matrix(0,ncol=2,nrow=length(genome[[locus]]$chrom1))
    for(i in 1:length(genome[[locus]]$chrom1)){
      all.config[i,] <- c(genome[[locus]]$chrom1[i],genome[[locus]]$chrom2[i])
    }
   return(all.config)
}

build.all.genotype <- function(genome){
  nb.locus <- get.nb.locus(genome)
  locus.all.config <- vector("list",nb.locus)
  for(locus in 1:nb.locus){
    locus.all.config[[locus]] <- build.genotype.from.locus(genome,locus)
  }
  all.haplotype <- build.all.haplotype(genome)
  all.genotype=c()
  nb.haplotype <- dim(all.haplotype)[1]
  for(haplotype1 in 1:nb.haplotype){
    for(haplotype2 in haplotype1:nb.haplotype){
      is.compatible <- TRUE
      for (locus in 1:nb.locus){
        allele.locus <- c(all.haplotype[haplotype1,locus],all.haplotype[haplotype2,locus])
        is.compatible <- is.compatible & genotype.is.in.locus(allele.locus,locus.all.config[[locus]])
      }
      if(is.compatible) all.genotype = c(all.genotype,haplotype1,haplotype2)
    }
  }
  all.genotype <- matrix(all.genotype,ncol=2,byrow = T)
  return(all.genotype)
}

#tell if a certain genotype configuration is in a locus
genotype.is.in.locus <- function(allele.locus,locus.all.config){
  is.here <- FALSE
  for (config in 1:dim(locus.all.config)[1]){
    is.here <- is.here | identical(allele.locus,locus.all.config[config,])
    is.here <- is.here | identical(allele.locus[2:1],locus.all.config[config,])
  }
  return(is.here)
}

#tells the position of a certain configuration
where.is.locus <- function(allele.locus,locus.all.config){
  where.is <- 0
  for (config in 1:dim(locus.all.config)[1]){
    if(all(allele.locus == locus.all.config[config,])) where.is <- config
    if(all(allele.locus[2:1] == locus.all.config[config,])) where.is <- config
  }
  if(where.is == 0) stop(paste("Cannot find",paste(allele.locus,collapse = " "), " in ",paste(locus.all.config,collapse = " "),"matdim:",paste(dim(locus.all.config),collapse = " ")))
  return(where.is)
}

get.genotype.index.from.haplotypes.index <- function(haplotypes,all.genotype){
  which(all.genotype[,1] == haplotypes[1] & all.genotype[,2]==haplotypes[2] | all.genotype[,1] == haplotypes[2] & all.genotype[,2]==haplotypes[1])
}

get.fitness.from.genotype.male <- function(genotype, genome,all.haplotype,all.genotype){
  fitness <- 1
  haplotype1 <- get.haplotype.from.index(all.genotype[genotype,1],get.nb.alleles.per.locus(genome))
  haplotype2 <- get.haplotype.from.index(all.genotype[genotype,2],get.nb.alleles.per.locus(genome))
  for (locus in 1:length(genome)){
    position <- where.is.locus(c(haplotype1[locus],haplotype2[locus]),build.genotype.from.locus(genome,locus))
    locus.fitness <- genome[[locus]]$fitness.male[position]
    fitness <- fitness*locus.fitness
  }
  return(fitness)
}

get.fitness.from.genotype.female <- function(genotype, genome,all.haplotype,all.genotype){
  fitness <- 1
  haplotype1 <- get.haplotype.from.index(all.genotype[genotype,1],get.nb.alleles.per.locus(genome))
  haplotype2 <- get.haplotype.from.index(all.genotype[genotype,2],get.nb.alleles.per.locus(genome))
  for (locus in 1:length(genome)){
    position <- where.is.locus(c(haplotype1[locus],haplotype2[locus]),build.genotype.from.locus(genome,locus))
    locus.fitness <- genome[[locus]]$fitness.female[position]
    fitness <- fitness*locus.fitness
  }
  return(fitness)
}
