
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
  nb.genotypes <- nrow(genome@all.genotype)
  nb.gamete <- nrow(genome@all.haplotype)
  males.genotype <- get.male(genome)
  maleness <- get.all.maleness(genome)
  male.frequency <- sum(initial.frequency[males.genotype]*maleness[males.genotype])
  relative.frequencies <- initial.frequency/male.frequency
  fitness.males <- get.all.fitness.male(genome)
  mean.fitness.male <- 0
  for(male.genotype in males.genotype) mean.fitness.male <- mean.fitness.male + 1/male.frequency*fitness.males[male.genotype]*maleness[male.genotype]*initial.frequency[male.genotype]

  relative.fitness.males <- fitness.males/mean.fitness.male
  gamete.frequency <- rep(0,nb.gamete)
  for(male.genotype in males.genotype){
    all.gamete <- get.gamete.and.frequency.from.genotype.male(genome,male.genotype)
    all.gamete$frequency <- all.gamete$frequency*
      relative.frequencies[male.genotype]*
      maleness[male.genotype]*
      relative.fitness.males[male.genotype]
    for(gamete in 1:length(all.gamete$frequency)){
      gamete.frequency[all.gamete$index[gamete]] <- gamete.frequency[all.gamete$index[gamete]] + all.gamete$frequency[gamete]
    }
  }
  return(gamete.frequency)
}


get.female.gamete.frequency <- function(genome,initial.frequency){
  nb.genotypes <- nrow(genome@all.genotype)
  nb.gamete <- nrow(genome@all.haplotype)
  females.genotype <- get.female(genome)
  femaleness <- get.all.femaleness(genome)
  female.frequency <- sum(initial.frequency[females.genotype]*femaleness[females.genotype])
  relative.frequencies <- initial.frequency/female.frequency
  fitness.females <- get.all.fitness.female(genome)
  mean.fitness.female <- 0
  for(female.genotype in females.genotype) mean.fitness.female <- mean.fitness.female + 1/female.frequency*fitness.females[female.genotype]*femaleness[female.genotype]*initial.frequency[female.genotype]

  relative.fitness.females <- fitness.females/mean.fitness.female
  gamete.frequency <- rep(0,nb.gamete)
  for(female.genotype in females.genotype){
    all.gamete <- get.gamete.and.frequency.from.genotype.female(genome,female.genotype)
    all.gamete$frequency <- all.gamete$frequency*
      relative.frequencies[female.genotype]*
      femaleness[female.genotype]*
      relative.fitness.females[female.genotype]
    for(gamete in 1:length(all.gamete$frequency)){
      gamete.frequency[all.gamete$index[gamete]] <- gamete.frequency[all.gamete$index[gamete]] + all.gamete$frequency[gamete]
    }
  }
  return(gamete.frequency)
}


