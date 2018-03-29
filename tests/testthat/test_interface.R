context("interface")

test_that("We can compute the evolution of the frequency through time",
          {
            freqs1 <- compute.frequency.evolution(genome1)
            freqs2 <- compute.frequency.evolution(genome.partially.sexual,generations=3)
            expect_error(compute.frequency.evolution(genome.partially.sexual,c(0.1,0.1,0.1),200))
            expect_equal(dim(freqs1),c(25,7))
            expect_equal(dim(freqs2),c(3,18))
            expect_equal(freqs1[3,],simulate.frequency(genome1,freqs1[2,]))
          }
)

