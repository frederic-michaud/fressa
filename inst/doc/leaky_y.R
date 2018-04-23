## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
locus1 = create.locus(allele1=c(1,1,1,2,3),
                      allele2 = c(1,2,3,3,3),
                      sd = c(0,1,0.9,1,1),
                      allele.name = c("x","y","L"))
locus2 = create.locus(allele1=  c(1,1,2),
                      allele2 = c(1,2,2),
                      fitness.female = c(1,0.9,0.8),
                      fitness.male = c(0.8,0.9,1),
                      allele.name = c("F","M"))
genome = create.genome(locus = list(locus1,locus2))
genome

## ------------------------------------------------------------------------
 #the third allele of the first locus appears with a frequency of 0.01
initial.frequency <- get.frequency.from.one.allele.frequency(genome,
                                                            locus = 1,
                                                             allele = 3,
                                                             allele.frequency =  0.01)
freqs <- compute.frequency.evolution(genome, initial.frequency, generations=1000)
plot.haplotype.frequency(genome,freqs)

## ------------------------------------------------------------------------
locus2 = create.locus(allele1=  c(1,1,2),
                      allele2 = c(1,2,2),
                      fitness.female = c(1,0.9,0.8),
                      fitness.male = c(0.6,0.8,1),
                      allele.name = c("F","M"))
genome = create.genome(locus = list(locus1,locus2))

## ------------------------------------------------------------------------

initial.frequency <- get.frequency.from.one.allele.frequency(genome,locus = 1,
                                                             allele = 3,
                                                             allele.frequency =  0.01)
freqs <- compute.frequency.evolution(genome, initial.frequency, generations=10000)
plot.haplotype.frequency(genome,freqs)

