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
locussd = create.locus(allele1=c(1,1),
                      allele2 = c(1,2),
                      sd = c(0,1),
                      fitness.male=c(1,1),
                      fitness.female=c(1,1),
                      allele.name = c("x","y"))
locusm = create.locus(allele1 = c(1,1,2),
                      allele2 = c(1,2,2),
                      recombination.modifier = c(0,0.5,1),
                      allele.name = c("-","+")
                      )
locussa = create.locus(allele1=  c(1,1,2),
                      allele2 = c(1,2,2),
                      fitness.female = c(0.6,1,0.6),
                      fitness.male = c(0.8,1,0.3),
                      allele.name = c("A","a"))
genome = create.genome(list(locussd,locussa,locusm),
                       male.recombination = c(0.01,0.01),
                       female.recombination = c(0.01,0.01))
print(genome)

## ------------------------------------------------------------------------
initial.frequency <- get.frequency.from.one.allele.frequency(genome,3,2,1.99)
freq <- compute.frequency.evolution(genome, generations = 1000)
plot.haplotype.frequency(genome,freq)

