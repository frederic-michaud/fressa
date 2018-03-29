context("frequency")

test_that("We get the correct kids from the parent (haplotype level)",
          {
            expect_equal(get.childs.genotype.from.parent.genotype(c(1,2), c(3,5)),matrix(c(1,1,2,2,3,5,3,5),ncol=2))
            expect_equal(get.childs.genotype.from.parent.genotype(c(1,1), c(3,5)),matrix(c(1,1,1,1,3,5,3,5),ncol=2))
          }
)

test_that("We get the correct kids from the parent (index level)",
          {
            expect_equal(get.childs.index.from.parent.index(2, 7,genome1),c(2,4,5,7))
            expect_equal(get.childs.index.from.parent.index(2, 7,genome2),c(1,7,2,24))
            expect_equal(get.childs.index.from.parent.index(2, 55,genome2),c(4,7,21,24))
          }
)

test_that("Frequencies are correctly updated at each generation",
          {
            expect_known_value(simulate.frequency(genome1,c(0.1,0.1,0.1,0.1,0.2,0.2,0.2)),"freq_genome1.rds",update = F)
            #expect_known_value(simulate.frequency(genome2,rep(1/171,171)),"freq_genome2.rds",update = F) #commented because too slow
            expect_known_value(simulate.frequency(genome.partially.sexual,rep(1/19,18)),"freq_genome_partially_sexual.rds",update = F)
          }
)
