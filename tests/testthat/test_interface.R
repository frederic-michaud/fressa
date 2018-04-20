context("interface")

test_that("We can generate the name of the haplotype",
          {
            haplotype.names <- get.haplotype.names(genome1)
            expect_equal(haplotype.names,c("xA","xa","yA","ya"))
            haplotype.names <- get.haplotype.names(genome2)
            expect_equal(haplotype.names,c("111","112","113","121","122","123","211","212","213","221","222","223","311","312","313","321","322","323"))
          }
)

test_that("We can generate the name of the genotype",
          {
            genotype.names <- get.genotype.names(genome1)
            expect_equal(genotype.names,c("xA|xA","xA|xa","xA|yA","xA|ya","xa|xa","xa|yA","xa|ya"))
          }
)

test_that("We can compute the evolution of the frequency through time",
          {
            freqs1 <- compute.frequency.evolution(genome1)
            freqs2 <- compute.frequency.evolution(genome.partially.sexual,generations=3)
            expect_error(compute.frequency.evolution(genome.partially.sexual,c(0.1,0.1,0.1),200))
            expect_equal(dim(freqs1),c(7,25))
            expect_equal(dim(freqs2),c(18,3))
            expect_equal(freqs1[,3],simulate.frequency(genome1,freqs1[,2]))
          }
)

test_that("We can compute haplotype frequency from genotype frequency",
          {
            freqs1 <- compute.frequency.evolution(genome1)
            expect_known_value(get.haplotype.frequency(genome1,freqs1),"freq_haplotype1.rds",update = F)
          }
)


test_that("We can compute allele frequency from genotype frequency",
          {
            freqs1 <- compute.frequency.evolution(genome1,generation=100)
            sd.allele <- get.allele.frequency(genome1,freqs1,1)
            sa.allele <- get.allele.frequency(genome1,freqs1,2)
            expect_equal(sd.allele[,100],c(0.75,0.25))
            expect_equal(sa.allele[,100],c(3/8-1/72,5/8+1/72),tolerance=6e-3)
          }
)



test_that("We can compute frequency of the genotypes from the frequency of one allele",
          {
            #The overall sum is one
            expect_equal(sum(get.frequency.from.one.allele.frequency(genome1,2,1,0.6)),1)
            expect_equal(sum(get.frequency.from.one.allele.frequency(genome2,1,2,0.3)),1)

            #the sum of the genotype which do have the allele is equal to what was given to the function
            expect_equal(sum(get.frequency.from.one.allele.frequency(genome1,2,1,0.1)[get.genotype.with.given.allele(genome1,2,1)]),0.1)
            expect_equal(sum(get.frequency.from.one.allele.frequency(genome2,1,1,0.3)[get.genotype.with.given.allele(genome2,1,1)]),0.3)
            expect_equal(sum(get.frequency.from.one.allele.frequency(genome.partially.sexual,2,1,0.71)[get.genotype.with.given.allele(genome.partially.sexual,2,1)]),0.71)
            expect_equal(sum(get.frequency.from.one.allele.frequency(genome.partially.sexual,2,1,0.25)[get.genotype.with.given.allele(genome.partially.sexual,2,1)]),0.25)
          }
)
