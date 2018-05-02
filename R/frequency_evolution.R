
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

get.male.gamete.frequency <- function(genome,initial.frequency){
  male.gamete.matrix <- genome@male.gamete.matrix
  gamete.frequency <- initial.frequency%*%male.gamete.matrix
  return(gamete.frequency/sum(gamete.frequency))
}


get.female.gamete.frequency <- function(genome,initial.frequency){
  female.gamete.matrix <- genome@female.gamete.matrix
  gamete.frequency <- initial.frequency%*%female.gamete.matrix
  return(gamete.frequency/sum(gamete.frequency))
}

build.male.gamete.matrix <- function(genome){
  nb.genotypes <- get.nb.genotype(genome)
  nb.gamete <- get.nb.gamete(genome)
  males.genotype <- get.male(genome)
  maleness <- get.all.maleness(genome)
  fitness.males <- get.all.fitness.male(genome)
  gamete.frequency <- rep(0,nb.gamete)
  list.of.gamete <- genome@all.gamete.male
  gamete.matrix <- as(matrix(0,nrow = nb.genotypes, ncol = nb.gamete),"sparseMatrix")
  for(male.genotype in males.genotype){
    all.gamete <- list.of.gamete[[male.genotype]]
    all.gamete$frequency <- all.gamete$frequency*
      maleness[male.genotype]*
      fitness.males[male.genotype]
    for(gamete in 1:length(all.gamete$frequency)){
      gamete.matrix[male.genotype, all.gamete$index[gamete]] <- gamete.matrix[male.genotype, all.gamete$index[gamete]] + all.gamete$frequency[gamete]
    }
  }
  return(gamete.matrix)
}


build.female.gamete.matrix <- function(genome){
  nb.genotypes <- get.nb.genotype(genome)
  nb.gamete <- get.nb.gamete(genome)
  females.genotype <- get.female(genome)
  femaleness <- get.all.femaleness(genome)
  fitness.females <- get.all.fitness.female(genome)
  gamete.frequency <- rep(0,nb.gamete)
  list.of.gamete <- genome@all.gamete.female
  gamete.matrix <- as(matrix(0,nrow = nb.genotypes, ncol = nb.gamete),"sparseMatrix")
  for(female.genotype in females.genotype){
    all.gamete <- list.of.gamete[[female.genotype]]
    all.gamete$frequency <- all.gamete$frequency*
      femaleness[female.genotype]*
      fitness.females[female.genotype]
    for(gamete in 1:length(all.gamete$frequency)){
      gamete.matrix[female.genotype, all.gamete$index[gamete]] <- gamete.matrix[female.genotype, all.gamete$index[gamete]] + all.gamete$frequency[gamete]
    }
  }
  return(gamete.matrix)
}
