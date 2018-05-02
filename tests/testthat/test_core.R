context("gamete")
test_that("number of locus is correctly retrieve from genome",
          {
            expect_equal(get.nb.locus(genome1),2)
            expect_equal(get.nb.locus(genome2),3)
          }
)

test_that("Number of genotype is correctly found",
          {
            expect_equal(get.nb.genotype(genome1),7)
            expect_equal(get.nb.genotype(genome2),18*19/2)
          }
)

test_that("An gamete is correctly computed from its index",
          {
            expect_equal(get.gamete.from.index(2,c(2,2)),c(1,2))
            expect_equal(get.gamete.from.index(16,c(3,3,2)),c(3,2,2))
          }
)

test_that("The number of allele per locus from a genome is correctly retreive",
          {
            expect_equal(get.nb.alleles.per.locus(genome1),c(2,2))
            expect_equal(get.nb.alleles.per.locus(genome2),c(3,2,3))
          }
)

test_that("All gamete are correctly found",
          {
            expect_equal(build.all.possible.gamete(genome1),matrix(c(1, 1, 2, 2, 1, 2, 1, 2),ncol=2))
            expect_equal(build.all.possible.gamete(genome2),matrix(c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3),ncol=3))
            }
)

test_that("All gamete are correctly found",
          {
            expect_equal(build.all.possible.gamete(genome1),matrix(c(1, 1, 2, 2, 1, 2, 1, 2),ncol=2))
            expect_equal(build.all.possible.gamete(genome2),matrix(c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3),ncol=3))
          }
)


test_that("We can build all genotype from a locus",
          {
            expect_equal(build.genotype.from.locus(genome1,1),matrix(c(1,1,1,2),ncol = 2))
            expect_equal(build.genotype.from.locus(genome1,2),matrix(c(1,1,2,1,2,2),ncol = 2))
            expect_equal(build.genotype.from.locus(genome2,1),matrix(c(1,1,1,2,2,3,1,2,3,2,3,3),ncol=2))
           }
)

test_that("We can build all possible genotype from a genome",
          {
            expect_equal(build.all.genotype(genome1),matrix(c(1, 1, 1, 1, 2, 2, 2, 1, 2, 3, 4, 2, 3, 4) ,ncol = 2))
            expect_equal(build.all.genotype(genome2),matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,13,13,13,13,13,13,14,14,14,14,14,15,15,15,15,16,16,16,17,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,5,6,7,8,9,10,11,12,13,14,15,16,17,18,6,7,8,9,10,11,12,13,14,15,16,17,18,7,8,9,10,11,12,13,14,15,16,17,18,8,9,10,11,12,13,14,15,16,17,18,9,10,11,12,13,14,15,16,17,18,10,11,12,13,14,15,16,17,18,11,12,13,14,15,16,17,18,12,13,14,15,16,17,18,13,14,15,16,17,18,14,15,16,17,18,15,16,17,18,16,17,18,17,18,18),ncol = 2))
          }
)

test_that("We can say if a locus configuration is in a set of locus configuration",
          {
            expect_true(genotype.is.in.locus(c(1,1),matrix(c(1,1,2,1,2,2),nrow = 3)))
            expect_true(genotype.is.in.locus(c(1,2),matrix(c(1,1,2,3,1,2,3,4),nrow = 4)))
            expect_true(genotype.is.in.locus(c(2,1),matrix(c(1,1,2,3,1,2,3,4),nrow = 4)))
            expect_false(genotype.is.in.locus(c(3,1),matrix(c(1,1,2,3,1,2,3,4),nrow = 4)))
          }
)

test_that("We can find the position of a locus configuration in a set of locus configuration ",
          {
            expect_equal(where.is.locus(c(1,1),matrix(c(1,1,2,1,2,2),nrow = 3)),1)
            expect_equal(where.is.locus(c(1,2),matrix(c(1,1,2,3,1,2,3,4),nrow = 4)),2)
            expect_equal(where.is.locus(c(2,1),matrix(c(1,1,2,3,1,2,3,4),nrow = 4)),2)
          }
)

test_that("Finding the index of a genotype from index of two gamete",
          {
            expect_equal(get.genotype.index.from.gametes.index(c(1,1),matrix(c(1,1,2,1,2,2),nrow = 3)),1)
            expect_equal(get.genotype.index.from.gametes.index(c(1,2),matrix(c(1,1,2,3,1,2,3,4),nrow = 4)),2)
            expect_equal(get.genotype.index.from.gametes.index(c(2,1),matrix(c(1,1,2,3,1,2,3,4),nrow = 4)),2)
          }
)

test_that("Fitness of male correctly computed",
          {
            expect_equal(get.fitness.from.genotype.male(2,genome1),0.8)
            expect_equal(get.fitness.from.genotype.male(3,genome1),0.6)
            expect_equal(get.fitness.from.genotype.male(2,genome2),0.66)
          }
)


test_that("Fitness of female correctly computed",
          {
            expect_equal(get.fitness.from.genotype.female(2,genome1),0.9)
            expect_equal(get.fitness.from.genotype.female(3,genome1),1)
            expect_equal(get.fitness.from.genotype.female(2,genome2),0.9)
          }
)

test_that("We can find the fitness of all males",
          {
            expect_equal(get.all.fitness.male(genome1),c(0.6, 0.8, 0.6, 0.8, 1.0, 0.8, 1.0))
            expect_equal(get.all.fitness.male(genome.partially.sexual),c(0.6,0.8,0.6,0.8,0.6,0.8,1,0.8,1,0.8,1,0.6,0.8,0.8,1,0.6,0.8,1))

          }
)

test_that("We can find the fitness of all females",
          {
            expect_equal(get.all.fitness.female(genome1),c(1, 0.9, 1, 0.9, 0.8, 0.9, 0.8))
            expect_equal(get.all.fitness.female(genome.partially.sexual),c(1,0.9,1,0.9,1,0.9,0.8,0.9,0.8,0.9,0.8,1,0.9,0.9,0.8,1,0.9,0.8))
          }
)
test_that("we can get all possible genotype which match a given gamete",{
    expect_equal(get.genotype.with.given.gamete(genome1,2),c(5,6,7,2,5))
    expect_equal(get.genotype.with.given.gamete(genome2,5),c(67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 5, 22, 38, 53, 67))
}
)

test_that("We can find all the gamete which contain a given allele",
          {
            expect_equal(get.gamete.with.given.allele(genome1,1,1),c(1,2))
            expect_equal(get.gamete.with.given.allele(genome1,2,2),c(2,4))
            expect_equal(get.gamete.with.given.allele(genome2,3,3),c(3,  6,  9, 12, 15, 18))
          }
)

test_that("We can find all the genotype which contain a given allele",
          {
            expect_equal(sort(get.genotype.with.given.allele(genome1,1,1)),c(1, 1, 2, 2, 3, 4, 5, 5, 6, 7))
            expect_equal(sort(get.genotype.with.given.allele(genome1,2,2)),c(2, 4, 5, 5, 6, 7, 7))
            expect_equal(sort(get.genotype.with.given.allele(genome2,1,3)),c(13,14,15,16,17,18,30,31,32,33,34,35,46,47,48,49,50,51,61,62,63,64,65,66,75,76,77,78,79,80,88,89,90,91,92,93,100,101,102,103,104,105,111,112,113,114,115,116,121,122,123,124,125,126,130,131,132,133,134,135,138,139,140,141,142,143,145,146,147,148,149,150,151,151,152,152,153,153,154,154,155,155,156,156,157,157,158,158,159,159,160,160,161,161,162,162,163,163,164,164,165,165,166,166,167,167,168,168,169,169,170,170,171,171))
          }
)

