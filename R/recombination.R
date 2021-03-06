
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

build.all.gamete <- function(genome,recombination.value){
  nb.genome <- get.nb.genotype(genome)
  recombination.list = as.list(1:nb.genome)
  recombination.modifier <- genome@all.recombination.modifier
  position.modifier <- genome@position.modifier
  for(genotype in 1:nb.genome)
  {
    scaled.recombination.value <- build.rescaled.recombination(recombination.value,
                                                               recombination.modifier,
                                                               genotype,
                                                               position.modifier)
    recombination.list[[genotype]] <- get.gamete.and.frequency.from.genotype.after.recombination(genome,
                                                                                                     genotype,
                                                                                                     scaled.recombination.value)
  }
  return(recombination.list)
}


get.gamete.and.frequency.from.genotype.after.recombination <- function(genome,genotype,recombination.value){
  all.genotype <- genome@all.genotype
  all.gamete <- genome@all.gamete
  if(length(recombination.value) == 0 ||sum(recombination.value) ==0){
    frequency <-  c(1/2,1/2)
    gamete.index <- all.genotype[genotype,]
  }
  else{
    gamete1.index <- all.genotype[genotype,1]
    gamete2.index <- all.genotype[genotype,2]
    gamete1 <- all.gamete[gamete1.index,]
    gamete2 <- all.gamete[gamete2.index,]
    nb.of.link <- length(recombination.value)
    frequency <- c()
    gamete.index <- c()
    for(recombination.index in 1:2^nb.of.link){
      gametes <- get.gamete.for.given.recombination(gamete1, gamete2, recombination.index)
      gamete.index.partial <-c(get.gamete.from.allele(genome,gametes[1,]),
                               get.gamete.from.allele(genome,gametes[2,])
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

get.gamete.for.given.recombination <- function(gamete1.before.recombination, gamete2.before.recombination, recombination.index){
  both.gamete = rbind(gamete1.before.recombination,gamete2.before.recombination)
  decomposition.recombination <- as.numeric(intToBits(recombination.index-1))
  nb.locus <- length(gamete1.before.recombination)
  gamete1 <- rep(0,nb.locus)
  gamete2 <- rep(0,nb.locus)
  gamete1 <- gamete1.before.recombination[1]
  gamete2 <- gamete2.before.recombination[1]
  for(i in 2:nb.locus){
    gamete1[i] <- both.gamete[mod(1+sum(decomposition.recombination[1:i-1]),2),i]
    gamete2[i] <- both.gamete[mod(sum(decomposition.recombination[1:i-1]),2),i]
  }
  return(rbind(gamete1,gamete2))
}

get.probability.for.given.recombination <- function(recombination.value, recombination.index){
  decomposition.recombination <- as.numeric(intToBits(recombination.index-1))
  nb.interval <- length(recombination.value)
  recombination.probability <- prod(abs(((1-recombination.value) - decomposition.recombination[1:nb.interval])))
  return(recombination.probability)
}

#return the recombination value after its scaling by the modifier
build.rescaled.recombination <- function(recombination.value,recombination.modifier,index.genome,position.modifier){
  scaled.recombination.value <- recombination.value
  if(length(recombination.modifier)>0)
    {
    if(length(position.modifier) == 0){
      scaled.recombination.value <- recombination.value*recombination.modifier[index.genome]
    }
    else{
      scaled.recombination.value[position.modifier] <- recombination.value[position.modifier]*recombination.modifier[index.genome]
    }
    }
  return(scaled.recombination.value)
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
