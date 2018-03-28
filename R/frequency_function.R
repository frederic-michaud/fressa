simulate.frequency <- function(genome,initial.frequency){
nb.genotypes <- dim(build.all.genotype(genome))[1]
males.genotype <- get.male(genome,build.all.haplotype(genome),build.all.genotype(genome))
females.genotype <- get.female(genome,build.all.haplotype(genome),build.all.genotype(genome))
nb.genotypes <- dim(build.all.genotype(genome))[1]
frequencies <- initial.frequency
new.frequencies <- rep(0,nb.genotypes)
male.frequency <- sum(frequencies[males.genotype])
female.frequency <- sum(frequencies[females.genotype])
props <- c()
all.genotype <- build.all.genotype(genome)
all.haplotype <- build.all.haplotype(genome)

fitness.males <- rep(0,nb.genotypes)
fitness.females <- rep(0,nb.genotypes)
for(i in 1:nb.genotypes) fitness.males[i] <- get.fitness.from.genotype.male(i,genome,all.haplotype, all.genotype)
for(i in 1:nb.genotypes) fitness.females[i] <- get.fitness.from.genotype.female(i,genome,all.haplotype, all.genotype)
mean.fitness.male <- 0
for(male.genotype in males.genotype) mean.fitness.male <- mean.fitness.male + 1/male.frequency*fitness.males[male.genotype]*frequencies[male.genotype]

mean.fitness.female <- 0
for(female.genotype in females.genotype) mean.fitness.female <- mean.fitness.female + 1/female.frequency*fitness.females[female.genotype]*frequencies[female.genotype]

  for(male.genotype in males.genotype){
    for(female.genotype in females.genotype){
      haplotype.male.1 <- all.genotype[male.genotype,1]
      haplotype.male.2 <- all.genotype[male.genotype,2]
      haplotype.female.1 <- all.genotype[female.genotype,1]
      haplotype.female.2 <- all.genotype[female.genotype,2]
      a1 = c(haplotype.male.1,haplotype.female.1)
      a2 = c(haplotype.male.1,haplotype.female.2)
      a3 = c(haplotype.male.2,haplotype.female.1)
      a4 = c(haplotype.male.2,haplotype.female.2)
      fitness.male <- fitness.males[male.genotype]/mean.fitness.male
      fitness.female <- fitness.females[female.genotype]/mean.fitness.female
      genotype.kid <- get.genotype.index.from.haplotypes.index(a1,all.genotype)
      new.frequencies[genotype.kid] <- new.frequencies[genotype.kid] +
        fitness.male*fitness.female*
        1/4*1/(male.frequency*female.frequency)*
        (frequencies[male.genotype]*frequencies[female.genotype])
      genotype.kid <- get.genotype.index.from.haplotypes.index(a2,all.genotype)
      new.frequencies[genotype.kid] <- new.frequencies[genotype.kid] +
        fitness.male*fitness.female*
        1/4*1/(male.frequency*female.frequency)*
        (frequencies[male.genotype]*frequencies[female.genotype])
      genotype.kid <- get.genotype.index.from.haplotypes.index(a3,all.genotype)
      new.frequencies[genotype.kid] <- new.frequencies[genotype.kid] +
        fitness.male*fitness.female*
        1/4*1/(male.frequency*female.frequency)*
        (frequencies[male.genotype]*frequencies[female.genotype])
      genotype.kid <- get.genotype.index.from.haplotypes.index(a4,all.genotype)
      new.frequencies[genotype.kid] <- new.frequencies[genotype.kid] +
        fitness.male*fitness.female*
        1/4*1/(male.frequency*female.frequency)*
        (frequencies[male.genotype]*frequencies[female.genotype])
    }
  }

new.frequencies
sum(new.frequencies)
male.frequency <- sum(new.frequencies[males.genotype])
female.frequency <- sum(new.frequencies[females.genotype])
male.frequency
female.frequency
frequencies <- new.frequencies
return(frequencies)
}



