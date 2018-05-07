context("interface subfunction")



test_that("We can compute correctly the marginal fitness of an allele",
          {
            init.freq <- get.frequency.from.one.allele.frequency(genome1,2,1,0.00000001)
            expect_equal(get.marginal.allele.fitness.single.generation.single.allele(genome1,init.freq,2,1),0.85)
            expect_equal(get.marginal.allele.fitness.single.generation.single.allele(genome1,init.freq,2,2),0.9)
          }
)

test_that("We can compute correctly the marginal fitness of all allele from one locus",
          {
            init.freq <- get.frequency.from.one.allele.frequency(genome1,2,1,0.00000001)
            expect_equal(get.marginal.allele.fitness.single.generation(genome1,init.freq,2),
                         c(get.marginal.allele.fitness.single.generation.single.allele(genome1,init.freq,2,1),
                           get.marginal.allele.fitness.single.generation.single.allele(genome1,init.freq,2,2)))
          }
)


test_that("We can correctly evalutate if convergence is reach",
          {
            freq1 <- matrix(c(-0.1,0.1,0.1,-0.05,0.05,0.05,-0.01,0.01,0.01),ncol=3)
            freq2 <- matrix(c(-0.1,-0.2,0.1,-0.05,-0.05,0.05,-0.01,-0.01,0.01),ncol=3)
            freq3 <- matrix(c(0.01,0.02,0.01,0.011,0.021,0.011,0.013,0.023,0.013),ncol=3)
            expect_true(is.converged(freq1, 3, 1, 1))
            expect_false(is.converged(freq1, 3, 1, 0.1))
            expect_false(is.converged(freq1, 2, 3, 1))
            expect_true(is.converged(freq2, 3, 1, 1))
            expect_false(is.converged(freq2, 3, 1, 0.1))
            expect_false(is.converged(freq2, 2, 3, 1))
            expect_false(is.converged(freq3, 3,1, 1))
          }
)
