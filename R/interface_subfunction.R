#' return names for the haplotype
get.haplotype.names <- function(genome){
  if(length(genome@locus[[1]]@allele.name) == 0) haplotype.names <- get.haplotype.names.from.allele.number(genome)
  else haplotype.names <- get.haplotype.names.from.allele.names(genome)
  return(haplotype.names)
}

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

plot.marginal.fitness <- function(genome,freqs,locus){
  allele.marginal.fitness <- get.marginal.allele.fitness(genome,freqs,locus)
  max.fit <- max(allele.marginal.fitness)
  min.fit <- min(allele.marginal.fitness)
  allele.number <- get.nb.alleles.per.locus(genome)[locus]
  palette <- get.palette(allele.number)
  plot(allele.marginal.fitness[1,],type="l",ylim=c(min.fit/1.2,1.2*max.fit),col=palette[1],xlab = "Generation",ylab="frequency")
  for(allele in 2:allele.number){
    lines(allele.marginal.fitness[allele,],col=palette[allele])
  }
  legend("topright",legend=get.allele.name(genome,locus),lty = rep(1,allele.number),col=palette)
}

get.marginal.allele.fitness <- function(genome,freqs,locus)
{
  nb.generation <- ncol(freqs)
  allele.number <- get.nb.alleles.per.locus(genome)[locus]
  allele.marginal.fitness <- matrix(0,ncol = nb.generation,nrow = allele.number)
  for (generation in 1:nb.generation){
    allele.marginal.fitness[,generation] <- get.marginal.allele.fitness.single.generation(genome,freqs[,generation],locus)
  }
  return(allele.marginal.fitness)
}


get.marginal.allele.fitness.single.generation <- function(genome,frequency,locus){
  nb.allele <- get.nb.alleles.per.locus(genome)[locus]
  marginal.fitnesses <- sapply(1:nb.allele,function(iter) get.marginal.allele.fitness.single.generation.single.allele(genome,frequency, locus,iter))
  return(marginal.fitnesses)
}

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

plot.haplotype.marginal.fitness <- function(genome,freqs){
  haplotype.marginal.fitness <- get.marginal.haplotype.fitness(genome,freqs)
  max.fit <- max(haplotype.marginal.fitness,na.rm = T)
  min.fit <- min(haplotype.marginal.fitness,na.rm = T)
  haplotype.number <- get.nb.haplotype(genome)
  palette <- get.palette(haplotype.number)
  plot(haplotype.marginal.fitness[1,],type="l",ylim=c(min.fit/1.2,1.2*max.fit),col=palette[1],xlab = "Generation",ylab="frequency")
  for(haplotype in 2:haplotype.number){
    lines(haplotype.marginal.fitness[haplotype,],col=palette[haplotype])
  }
  legend("topright",legend=get.haplotype.names(genome),lty = rep(1,haplotype.number),col=palette)
}

get.marginal.haplotype.fitness <- function(genome,freqs)
{
  nb.generation <- ncol(freqs)
  haplotype.number <- get.nb.haplotype(genome)
  haplotype.marginal.fitness <- matrix(0,ncol = nb.generation,nrow = haplotype.number)
  for (generation in 1:nb.generation){
    haplotype.marginal.fitness[,generation] <- get.marginal.haplotype.fitness.single.generation(genome,freqs[,generation])
  }
  return(haplotype.marginal.fitness)
}


get.marginal.haplotype.fitness.single.generation <- function(genome,frequency){
  nb.haplotype <- get.nb.haplotype(genome)
  marginal.fitnesses <- sapply(1:nb.haplotype,function(iter) get.marginal.haplotype.fitness.single.generation.single.haplotype(genome,frequency, iter))
  return(marginal.fitnesses)
}

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
