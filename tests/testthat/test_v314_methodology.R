test_that("v3.14 TME adjustment grid includes the over-adjusted negative control", {
  specs <- tme_adjustment_specs()
  expect_true(all(c("S0", "S1", "S2", "S3", "S4") %in% names(specs)))
  expect_equal(specs$S1, c("stromal_caf", "immune_pan"))
  expect_true("ductal_epithelial" %in% specs$S4)
  expect_true("acinar_pancreas" %in% specs$S4)
})

test_that("v3.14 identifier audit flags noncoding or legacy transition symbols", {
  tbl <- tibble::tibble(gene = c("TP53", "LOC100130175", "UCKL1-AS1", "CRYBB2P1"))
  out <- annotate_transition_identifier_class(tbl)
  expect_true("transition_external_validation_tier" %in% names(out))
  expect_true(any(out$noncoding_or_legacy_symbol_like))
})

test_that("v3.14 ordinary limma SE extraction is independent of moderated t", {
  skip_if_not_installed("limma")
  set.seed(1)
  expr <- matrix(rnorm(40), nrow = 10)
  rownames(expr) <- paste0("G", seq_len(nrow(expr)))
  colnames(expr) <- paste0("S", seq_len(ncol(expr)))
  design <- stats::model.matrix(~ factor(c(0, 0, 1, 1)))
  fit <- limma::lmFit(expr, design)
  se <- ordinary_se_from_limma(fit, 2)
  expect_equal(length(se), nrow(expr))
  expect_true(all(is.finite(se)))
})
