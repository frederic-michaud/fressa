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
locus2 = create.locus(allele1=  c(1,1,2),
                      allele2 = c(1,2,2),
                      fitness.female = c(1,0.9,0.8),
                      fitness.male = c(0.8,0.9,1),
                      allele.name = c("F","M"))
genome = create.genome(list(locus1,locus2))
print(genome)

## ------------------------------------------------------------------------
freqs <- compute.frequency.evolution(genome, generations=200)

## ------------------------------------------------------------------------
plot.haplotype.frequency(genome, freqs)

## ------------------------------------------------------------------------
#locus 1 is the same as before
locus2 = create.locus(allele1=  c(1,1,2),
                      allele2 = c(1,2,2),
                      fitness.female = c(1,0.9,0.8),
                      fitness.male = c(0.6,0.8,1),
                      allele.name = c("F","M"))
genome = create.genome(list(locus1,locus2))
print(genome)

## ------------------------------------------------------------------------
freqs <- compute.frequency.evolution(genome, generations=200)
plot.haplotype.frequency(genome, freqs)

## ------------------------------------------------------------------------
all.final.polymorphism <- c()
sequence <- seq(0,2,0.025)
for(assymetry in sequence)
{
locus2 = create.locus(allele1=  c(1,1,2),
                      allele2 = c(1,2,2),
                      fitness.female = c(1,0.9,0.8),
                      fitness.male = c(0.8-2*assymetry*0.1,0.9-assymetry*0.1,1),
                      allele.name = c("F","M"))
genome = create.genome(list(locus1,locus2))
freqs <- compute.frequency.evolution(genome, generations=500)
haplotype.freq <- get.haplotype.frequency(genome, freqs)
final.polymorphism <- haplotype.freq[1,ncol(haplotype.freq)]-haplotype.freq[2,ncol(haplotype.freq)]
all.final.polymorphism <- c(all.final.polymorphism,final.polymorphism)
}

## ------------------------------------------------------------------------
plot(sequence,all.final.polymorphism,xlab = "strength of assymetry",ylab="xm-xf",main="evolution of the polymorphism when changing the SA")

## ------------------------------------------------------------------------
locus1 = create.locus(allele1=c(1,1),
                      allele2 = c(1,2),
                      sd = c(0,1),
                      fitness.male=c(1,1),
                      fitness.female=c(1,1),
                      allele.name = c("x","y"))
locus2 = create.locus(allele1=  c(1,1,2),
                      allele2 = c(1,2,2),
                      fitness.female = c(1,0.9,0.8),
                      fitness.male = c(0.8,0.9,1),
                      allele.name = c("F","M"))
genome = create.genome(list(locus1,locus2),
                       male.recombination = c(0.05),
                       female.recombination = c(0.05))

## ------------------------------------------------------------------------
freqs <- compute.frequency.evolution(genome, generations=200)
plot.haplotype.frequency(genome, freqs)

## ------------------------------------------------------------------------
all.final.polymorphism <- c()
all.recombination.rate <- seq(0,0.5,0.01)
for(recombination.rate in all.recombination.rate){
  genome = create.genome(list(locus1,locus2),male.recombination = c(recombination.rate),female.recombination = c(recombination.rate))
  freqs <- compute.frequency.evolution(genome, generations=200)
  haplotype.freq <- get.haplotype.frequency(genome, freqs)
  final.polymorphism <- haplotype.freq[1,ncol(haplotype.freq)]-haplotype.freq[2,ncol(haplotype.freq)]
  all.final.polymorphism <- c(all.final.polymorphism,final.polymorphism)
}
plot(all.recombination.rate,all.final.polymorphism,xlab = "recombination rate",ylab = "Xm-Xf","polymorphism as a function of the recombination rate")

