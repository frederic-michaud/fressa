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
  all.genotype <- build.all.genotype(genome)
  childs <- get.childs.genotype.from.parent.genotype(all.genotype[male.genotype,],all.genotype[female.genotype,])
  a1 <- get.genotype.index.from.haplotypes.index(childs[1,],all.genotype)
  a2 <- get.genotype.index.from.haplotypes.index(childs[2,],all.genotype)
  a3 <- get.genotype.index.from.haplotypes.index(childs[3,],all.genotype)
  a4 <- get.genotype.index.from.haplotypes.index(childs[4,],all.genotype)
  return(c(a1,a2,a3,a4))
}

#' Get the frequency of new born as a function of adults

simulate.frequency <- function(genome,initial.frequency){
nb.genotypes <- dim(build.all.genotype(genome))[1]
males.genotype <- get.male(genome)
females.genotype <- get.female(genome)
frequencies <- initial.frequency
male.frequency <- sum(frequencies[males.genotype])
female.frequency <- sum(frequencies[females.genotype])
props <- c()
all.genotype <- build.all.genotype(genome)
all.haplotype <- build.all.haplotype(genome)
new.frequencies <- rep(0,nb.genotypes)
fitness.males <- sapply(1:nb.genotypes, get.fitness.from.genotype.male,genome = genome)
fitness.females <- sapply(1:nb.genotypes, get.fitness.from.genotype.female,genome = genome)

mean.fitness.male <- 0
for(male.genotype in males.genotype) mean.fitness.male <- mean.fitness.male + 1/male.frequency*fitness.males[male.genotype]*frequencies[male.genotype]
fitness.males <- fitness.males/mean.fitness.male
mean.fitness.female <- 0
for(female.genotype in females.genotype) mean.fitness.female <- mean.fitness.female + 1/female.frequency*fitness.females[female.genotype]*frequencies[female.genotype]
fitness.females <- fitness.females/mean.fitness.female

for(male.genotype in males.genotype){
  for(female.genotype in females.genotype){
    childs <- get.childs.index.from.parent.index(male.genotype,female.genotype,genome)
    fitness.male <- fitness.males[male.genotype]
    fitness.female <- fitness.females[female.genotype]
    fitness.couple <- fitness.male*fitness.female
    probability.of.couple <-1/(male.frequency*female.frequency)*(frequencies[male.genotype]*frequencies[female.genotype])
    for(child in childs){
      new.frequencies[child] <- new.frequencies[child] + fitness.couple*probability.of.couple/4
    }
  }
}
return(new.frequencies)
}



