context("interface")

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

test_that("We can find all the haplotype which contain a given allele",
          {
            expect_equal(get.haplotype.with.given.allele(genome1,1,1),c(1,2))
            expect_equal(get.haplotype.with.given.allele(genome1,2,2),c(2,4))
            expect_equal(get.haplotype.with.given.allele(genome2,3,3),c(3,  6,  9, 12, 15, 18))
          }
)

test_that("We can find all the genotype which contain a given allele",
          {
            expect_equal(sort(get.genotype.with.given.allele(genome1,1,1)),c(1, 1, 2, 2, 3, 4, 5, 5, 6, 7))
            expect_equal(sort(get.genotype.with.given.allele(genome1,2,2)),c(2, 4, 5, 5, 6, 7, 7))
            expect_equal(sort(get.genotype.with.given.allele(genome2,1,3)),c(13,14,15,16,17,18,30,31,32,33,34,35,46,47,48,49,50,51,61,62,63,64,65,66,75,76,77,78,79,80,88,89,90,91,92,93,100,101,102,103,104,105,111,112,113,114,115,116,121,122,123,124,125,126,130,131,132,133,134,135,138,139,140,141,142,143,145,146,147,148,149,150,151,151,152,152,153,153,154,154,155,155,156,156,157,157,158,158,159,159,160,160,161,161,162,162,163,163,164,164,165,165,166,166,167,167,168,168,169,169,170,170,171,171))
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
