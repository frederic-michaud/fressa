context("modifier")

test_that("We can correctly build the recombination for any genotype",
          {
            expect_equal(get.all.recombination.modifier(genome.with.modifier),c(1,1,0,0,0,0,1,0,0,0,0,1,1,1,1,1,1,1))
          }
)
