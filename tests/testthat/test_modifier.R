context("modifier")

test_that("We can correctly build the recombination for any genotype",
          {
            expect_equal(get.all.recombination.modifier(genome.with.modifier),c(1,1,0,0,0,0,1,0,0,0,0,1,1,1,1,1,1,1))
          }
)


build.rescaled.recombination

test_that("We can compute a new recombination rate taking into account modifier",
                                      {
                                        #multiply by 1 (third value(= index.genome) of recombination.modifier) the second (position.modifier)
                                        expect_equal(build.rescaled.recombination(recombination.value = c(0.1,0.1,0.1),
                                                                                  recombination.modifier = c(1,2,1,2), #four genotype
                                                                                  index.genome = 3,
                                                                                  position.modifier = c(2)
                                                                                  ),c(0.1,0.1,0.1))
                                        expect_equal(build.rescaled.recombination(recombination.value = c(0.1,0.1,0.1),
                                                                                  recombination.modifier = c(1,2,2,2),
                                                                                  index.genome = 3,
                                                                                  position.modifier = c(2)
                                        ),c(0.1,0.2,0.1))
                                        expect_equal(build.rescaled.recombination(recombination.value = c(0.15,0.1,0.1),
                                                                                  recombination.modifier = c(3,2,3,2),
                                                                                  index.genome = 1,
                                                                                  position.modifier = c()
                                        ),c(0.45,0.3,0.3))
                                        expect_equal(build.rescaled.recombination(recombination.value = c(0.1,0.1,0.1),
                                                                                  recombination.modifier = c(3,1.5,3,2),
                                                                                  index.genome = 2,
                                                                                  position.modifier = c(1,2)
                                        ),c(0.15,0.15,0.1))
                                      }
)
