test_that("clean_symbol handles multi-mapping", {
  expect_equal(clean_symbol("TP53 /// ABC"), "TP53")
  expect_true(is.na(clean_symbol("---")))
})

test_that("signature definitions are non-empty", {
  sigs <- builtin_tme_signatures()
  expect_true(length(sigs) >= 4)
  expect_true(all(lengths(sigs) >= 5))
})
