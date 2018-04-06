context("interface subfunction")

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
