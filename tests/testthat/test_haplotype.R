context("haplotype")
test_that("number of haplotype is correctly computed",
          {
            expect_equal(get.nb.genotype(genome1),7)
            expect_equal(get.nb.genotype(genome2),18*19/2)
          }
)

