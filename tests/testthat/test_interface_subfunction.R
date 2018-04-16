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
