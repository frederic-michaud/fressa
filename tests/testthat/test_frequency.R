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
