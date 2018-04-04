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
