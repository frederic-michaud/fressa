context("sexual")

test_that("Id of sd is correctly found",
          {
            expect_equal(get.id.sd.locus(genome1),1)
            expect_equal(get.id.sd.locus(genome2),3)
          }
)

test_that("If the number of sd is different from one, get error",
          {
            expect_error(get.id.sd.locus(create.genome(locus = list(locus2,locus2))))
            expect_error(get.id.sd.locus(create.genome(locus = list(locus3,locus3))))
          }
)

test_that("The function is.male correctly identified male",
          {
            expect_true(is.male(3,genome1))
            expect_false(is.male(5,genome1))
          }
)

test_that("The function is.female correctly identified female",
          {
            expect_false(is.female(3,genome1))
            expect_true(is.female(5,genome1))
          }
)

test_that("The function get.male get all the male present in the population",
          {
            expect_equal(get.male(genome1),c(3,4,6,7))
            expect_equal(get.male(genome2),c(2,3,5,6,8,9,11,12,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,53,54,56,57,59,60,62,63,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,95,96,98,99,101,102,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,128,129,131,132,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,152,153,155,156,157,158,159,160,161,162,163,164,165,167,168,169,170,171))
            expect_equal(get.male(genome.partially.sexual),c(3,4,8,9,12,13,14,15,16,17,18))
          }
)

test_that("The function get.all.maleness getthe correct value fot the maleness",
          {
            expect_equal(get.all.maleness(genome1),c(0,0,1,1,0,1,1))
            expect_equal(get.all.maleness(genome.partially.sexual),c(0,0,1,1,0,0,0,1,1,0,0,1,1,1,1,0.5,0.5,0.5))
          }
)
test_that("The function get.maleness correctly get the maleness",
          {
            expect_equal(get.maleness(3,genome1),1)
            expect_equal(get.maleness(5,genome1),0)
            expect_equal(get.maleness(17,genome.partially.sexual),0.5)
          }
)

test_that("The function get.all.femaleness getthe correct value fot the femaleness",
          {
            expect_equal(get.all.femaleness(genome1),1-c(0,0,1,1,0,1,1))
            expect_equal(get.all.femaleness(genome.partially.sexual),1-c(0,0,1,1,0,0,0,1,1,0,0,1,1,1,1,0.5,0.5,0.5))
          }
)
test_that("The function get.femaleness correctly get the femaleness",
          {
            expect_equal(get.femaleness(3,genome1),0)
            expect_equal(get.femaleness(5,genome1),1)
            expect_equal(get.femaleness(17,genome.partially.sexual),0.5)

          }
)

test_that("The function get.gamete.male & get.gamete.female correctly return the list of all gamete present in male, resp. female",
          {
            expect_equal(get.gamete.male(genome1),c(1,2,3,4))
            expect_equal(get.gamete.female(genome1),c(1,2))
            expect_equal(get.gamete.male(genome.partially.sexual),1:6)
            expect_equal(get.gamete.female(genome.partially.sexual),c(1,2,5,6))
          }
)
