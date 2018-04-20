context("Utilities")
test_that("test the modulo function whith n%n = n instead of 0",
          {
            expect_equal(mod(5,3),2)
            expect_equal(mod(6,3),3)
            expect_equal(mod(37,3),1)
            expect_equal(mod(-2,3),1)
            expect_equal(mod(11,13),11)
            expect_equal(mod(13,13),13)
            expect_equal(mod(19,13),6)
          }
)

