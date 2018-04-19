context("frequency")

test_that("Frequencies are correctly updated at each generation",
          {
            expect_known_value(simulate.frequency(genome1,c(0.1,0.1,0.1,0.1,0.2,0.2,0.2)),"freq_genome1.rds",update = F)
            expect_known_value(simulate.frequency(genome2,rep(1/171,171)),"freq_genome2.rds",update = F) #can be commented if too slow.
            expect_known_value(simulate.frequency(genome.partially.sexual,rep(1/19,18)),"freq_genome_partially_sexual.rds",update = F)
          }
)
