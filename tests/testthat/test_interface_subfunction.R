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
