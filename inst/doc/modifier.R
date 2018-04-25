## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
locus1 = create.locus(allele1=c(1,1),
                      allele2 = c(1,2),
                      sd = c(0,1),
                      fitness.male=c(1,1),
                      fitness.female=c(1,1),
                      allele.name = c("x","y"))
locus2 = create.locus(allele1 = c(1,1,2),
                      allele2 = c(1,2,2),
                      recombination.modifier = c(0,0.5,1),
                      allele.name = c("-","+")
                      )
locus3 = create.locus(allele1=  c(1,1,2),
                      allele2 = c(1,2,2),
                      fitness.female = c(0.9,0.95,1),
                      fitness.male = c(1,0.95,0.9),
                      allele.name = c("M","F"))
genome = create.genome(list(locus1,locus2,locus3),
                       male.recombination = c(0,0.01),
                       female.recombination = c(0,0.01))
print(genome)

## ------------------------------------------------------------------------
freq <- compute.frequency.evolution(genome,generations = 2500)
plot.haplotype.frequency(genome,freq)

## ------------------------------------------------------------------------
locussa = create.locus(allele1=  c(1,1,2),
                      allele2 = c(1,2,2),
                      fitness.female = c(0.75,1,0.75),
                      fitness.male = c(0.95,1,0.4),
                      allele.name = c("a","A"))
genome = create.genome(list(locussd,locussa,locusm),
                       male.recombination = c(0.01,0.5),
                       female.recombination = c(0.01,0.5))
print(genome)
get.haplotype.names(genome)

## ------------------------------------------------------------------------
initial.frequency <- get.frequency.from.gamete.frequency(genome,
                                                 male.gamete.frequency = c(0.49,0,0.01,0,0.01,0,0.49,0),
                                                 female.gamete.frequency = c(0.99,0,0.01,0,0,0,0,0))
freq <- compute.frequency.evolution(genome,initial.frequency, generations = 1000)
plot.haplotype.frequency(genome,freq)
equilibrium.frequency <- freq[,ncol(freq)]
equilibrium.frequency[8] <- 0.01 #addinf a mutation

## ------------------------------------------------------------------------
freq <- compute.frequency.evolution(genome,equilibrium.frequency, generations = 40000)
plot.haplotype.frequency(genome,freq)
plot.allele.frequency(genome,freq,3)

## ------------------------------------------------------------------------
locussa = create.locus(allele1=  c(1,1,2),
                      allele2 = c(1,2,2),
                      fitness.female = c(0.75,1,0.75),
                      fitness.male = c(0.9,1,0.6),
                      allele.name = c("a","A"))
genome = create.genome(list(locussd,locussa,locusm),
                       male.recombination = c(0.01,0.5),
                       female.recombination = c(0.01,0.5))
print(genome)
get.genotype.names(genome)

## ------------------------------------------------------------------------
initial.frequency <- get.frequency.from.gamete.frequency(genome,
                                                 male.gamete.frequency = c(0.49,0,0.01,0,0.01,0,0.49,0),
                                                 female.gamete.frequency = c(0.99,0,0.01,0,0,0,0,0))
#initial.frequency <- get.frequency.from.gamete.frequency(genome,
#                                                male.gamete.frequency = c(0.49,0,0.01,0,0.49,0,0.01,0),
#                                               female.gamete.frequency = c(0.99,0,0.01,0,0,0,0,0))
freq <- compute.frequency.evolution(genome,initial.frequency, generations = 1000)
plot.haplotype.frequency(genome,freq)
equilibrium.frequency <- freq[,ncol(freq)]
equilibrium.frequency[8] <- 0.01 #adding a mutation

## ------------------------------------------------------------------------
freq <- compute.frequency.evolution(genome,equilibrium.frequency, generations = 40000)
plot.haplotype.frequency(genome,freq)
plot.allele.frequency(genome,freq,3)

## ------------------------------------------------------------------------
locussa = create.locus(allele1=  c(1,1,2),
                      allele2 = c(1,2,2),
                      fitness.female = c(0.2,1,0.2),
                      fitness.male = c(0.2,1,0.2),
                      allele.name = c("a","A"))
genome = create.genome(list(locussd,locussa,locusm),
                       male.recombination = c(0.01,0.5),
                       female.recombination = c(0.01,0.5))
print(genome)
get.genotype.names(genome)

## ------------------------------------------------------------------------
initial.frequency <- get.frequency.from.gamete.frequency(genome,
                                                 male.gamete.frequency = c(0.49,0,0.01,0,0.01,0,0.49,0),
                                                 female.gamete.frequency = c(0.99,0,0.01,0,0,0,0,0))
#initial.frequency <- get.frequency.from.gamete.frequency(genome,
#                                                male.gamete.frequency = c(0.49,0,0.01,0,0.49,0,0.01,0),
#                                               female.gamete.frequency = c(0.99,0,0.01,0,0,0,0,0))
freq <- compute.frequency.evolution(genome,initial.frequency, generations = 1000)
plot.haplotype.frequency(genome,freq)
equilibrium.frequency <- freq[,ncol(freq)]
equilibrium.frequency[8] <- 0.025 #adding a mutation

## ------------------------------------------------------------------------
freq <- compute.frequency.evolution(genome,equilibrium.frequency, generations = 10000)
plot.haplotype.frequency(genome,freq)
plot.allele.frequency(genome,freq,3)

