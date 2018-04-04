context("interface subfunction")

test_that("We can generate the name of the haplotype",
          {
            haplotype.names <- get.haplotype.names(genome1)
            expect_equal(haplotype.names,c("11","12","21","22"))
          }
)

test_that("We can generate the name of the genotype",
          {
            genotype.names <- get.genotype.names(genome1)
            expect_equal(genotype.names,c("11|11", "11|12", "11|21", "11|22", "12|12", "12|21", "12|22"))
          }
)
