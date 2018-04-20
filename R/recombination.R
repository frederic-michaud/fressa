
get.gamete.and.frequency.from.genotype.male <- function(genome,genotype){
  recombination.value <- genome@male.recombination
  male.gamete <- get.gamete.and.frequency.from.genotype.after.recombination(genome,genotype,recombination.value)
  return(male.gamete)
}

get.gamete.and.frequency.from.genotype.female <- function(genome,genotype){
  recombination.value <- genome@female.recombination
  female.gamete <- get.gamete.and.frequency.from.genotype.after.recombination(genome,genotype,recombination.value)
  return(female.gamete)
}

get.gamete.and.frequency.from.genotype.after.recombination <- function(genome,genotype,recombination.value){
  all.genotype <- genome@all.genotype
  all.haplotype <- genome@all.haplotype
  if(length(recombination.value) ==0){
    frequency <-  c(1/2,1/2)
    gamete.index <- all.genotype[genotype,]
  }
  else{
    haplotype1.index <- all.genotype[genotype,1]
    haplotype2.index <- all.genotype[genotype,2]
    haplotype1 <- all.haplotype[haplotype1.index,]
    haplotype2 <- all.haplotype[haplotype2.index,]
    nb.of.link <- length(recombination.value)
    frequency <- c()
    gamete.index <- c()
    for(recombination.index in 1:2^nb.of.link){
      gametes <- get.gamete.for.given.recombination(haplotype1, haplotype2, recombination.index)
      gamete.index.partial <-c(get.haplotype.from.allele(genome,gametes[1,]),
                               get.haplotype.from.allele(genome,gametes[2,])
      )
      gamete.frequency.partial <- get.probability.for.given.recombination(recombination.value, recombination.index)
      #gamete.frequency.partial is twice in next expression, once for each gamete.
      frequency <- c(frequency,gamete.frequency.partial/2,gamete.frequency.partial/2)
      gamete.index <- c(gamete.index,gamete.index.partial)
    }
  }
  gamete.with.frequency <- data.frame(frequency = frequency,index = gamete.index)
  gamete.with.frequency <- glue.frequency(gamete.with.frequency)
  return(gamete.with.frequency)
}

get.gamete.for.given.recombination <- function(haplotype1, haplotype2, recombination.index){
  both.haplotype = rbind(haplotype1,haplotype2)
  decomposition.recombination <- as.numeric(intToBits(recombination.index-1))
  nb.locus <- length(haplotype1)
  gamete1 <- rep(0,nb.locus)
  gamete2 <- rep(0,nb.locus)
  gamete1 <- haplotype1[1]
  gamete2 <- haplotype2[1]
  for(i in 2:nb.locus){
    gamete1[i] <- both.haplotype[mod(1+sum(decomposition.recombination[1:i-1]),2),i]
    gamete2[i] <- both.haplotype[mod(sum(decomposition.recombination[1:i-1]),2),i]
  }
  return(rbind(gamete1,gamete2))
}

get.probability.for.given.recombination <- function(recombination.value, recombination.index){
  decomposition.recombination <- as.numeric(intToBits(recombination.index-1))
  nb.interval <- length(recombination.value)
  recombination.probability <- prod(abs(((1-recombination.value) - decomposition.recombination[1:nb.interval])))
  return(recombination.probability)
}

build.all.gamete <- function(genome,recombination.value){
nb.genome <- get.nb.genotype(genome)
  recombination.list = as.list(1:nb.genome)
  for(i in 1:nb.genome)
  {
    recombination.list[[i]] <- get.gamete.and.frequency.from.genotype.after.recombination(genome,i,recombination.value)
  }
  return(recombination.list)
}

glue.frequency <- function(frequency){
  nb.gamete <- length(frequency$index)
  max.gamete <- max(frequency$index+1)
  frequency.vector <- rep(0,max.gamete)
  for(i in 1:nb.gamete){
    frequency.vector[frequency$index[i]] <- frequency.vector[frequency$index[i]] + frequency$frequency[i]
  }
  new.frequency <- c()
  new.index <- c()
  for(i in 1:max.gamete){
    if(frequency.vector[i] > 0){
      new.frequency <- c(new.frequency,frequency.vector[i])
      new.index <- c(new.index,i)
    }
  }
return(data.frame(frequency = new.frequency,index = new.index))
}
