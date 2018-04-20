context("Recombination")

test_that("we can get all possible gamete from a male",{
  gamete <- get.gamete.and.frequency.from.genotype.male(genome1,2)
  expect_equal(gamete$frequency,c(0.5,0.5))
  expect_equal(gamete$index,c(1,2))
  gamete <- get.gamete.and.frequency.from.genotype.male(genome2,7)
  expect_equal(gamete$frequency,c(0.5,0.5))
  expect_equal(gamete$index,c(1,7))
  gamete <- get.gamete.and.frequency.from.genotype.male(genome.with.recomb,7)
  expect_equal(gamete$frequency,c(0.5,0.5))
  expect_equal(gamete$index,c(2,4))
  gamete <- get.gamete.and.frequency.from.genotype.male(genome.with.recomb,4)
  expect_equal(gamete$frequency,c(0.45,0.05,0.05,0.45))
  expect_equal(gamete$index,c(1,2,3,4))
}
)


test_that("we can get all possible gamete from a female",{
  gamete <- get.gamete.and.frequency.from.genotype.female(genome1,2)
  expect_equal(gamete$frequency,c(0.5,0.5))
  expect_equal(gamete$index,c(1,2))
  gamete <- get.gamete.and.frequency.from.genotype.female(genome2,7)
  gamete <- get.gamete.and.frequency.from.genotype.female(genome.with.recomb,4)
  expect_equal(gamete$frequency,c(0.475,0.025,0.025,0.475))
  expect_equal(gamete$index,c(1,2,3,4))
}
)

test_that("We can compute the gamete from two haplotype and a recombination index",{
  expect_equal(get.gamete.for.given.recombination(c(1,2,3,4),c(5,6,7,8),1),rbind(gamete1 = c(1,2,3,4),gamete2 = c(5,6,7,8)))
  expect_equal(get.gamete.for.given.recombination(c(1,2,3,4),c(5,6,7,8),2),rbind(gamete1 = c(1,6,7,8),gamete2 = c(5,2,3,4)))
  expect_equal(get.gamete.for.given.recombination(c(1,2,3,4),c(5,6,7,8),3),rbind(gamete1 = c(1,2,7,8),gamete2 = c(5,6,3,4)))
  expect_equal(get.gamete.for.given.recombination(c(1,2,3,4),c(5,6,7,8),4),rbind(gamete1 = c(1,6,3,4),gamete2 = c(5,2,7,8)))
}
)

test_that("We can compute which gamete appears with which probability in recombination",{
  expect_equal(get.probability.for.given.recombination(c(0.1,0.2,0.3),1),0.9*0.8*0.7)
  expect_equal(get.probability.for.given.recombination(c(0.1,0.2,0.3),2),0.1*0.8*0.7)
  expect_equal(get.probability.for.given.recombination(c(0.1,0.2,0.3),3),0.9*0.2*0.7)
  expect_equal(get.probability.for.given.recombination(c(0.1,0.2,0.3),4),0.1*0.2*0.7)
  expect_equal(get.probability.for.given.recombination(c(0.1,0.2,0.3),8),0.1*0.2*0.3)
}
)


test_that("we can get build correctly the list of all possible gamete",{
  expect_known_value(build.all.gamete(genome.with.recomb,genome.with.recomb@female.recombination),"recomb_genome_female.rds",update = F)
  expect_known_value(build.all.gamete(genome.with.recomb,genome.with.recomb@male.recombination),"recomb_genome_male.rds",update = F)
  #this was checked by hand before being saved!
  expect_known_value(build.all.gamete(genome.with.modifier,genome.with.modifier@male.recombination),"recomb_modifier_genome_male.rds",update = F)
}
)
