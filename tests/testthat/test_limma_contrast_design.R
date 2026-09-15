test_that("benign_vs_normal uses control as reference and benign as coefficient", {
  meta <- tibble::tibble(
    sample = c("S1", "S2", "S3", "S4"),
    condition = c("control", "control", "benign", "benign")
  )
  cf <- condition_factor_for_contrast(meta, "benign_vs_normal")
  m <- meta[cf$keep, , drop = FALSE]
  m$condition_raw <- m$condition
  condition <- droplevels(factor(
    ifelse(m$condition_raw == cf$case_level, cf$case_level, cf$reference_level),
    levels = c(cf$reference_level, cf$case_level)
  ))
  m$condition <- condition
  design <- stats::model.matrix(~ condition, data = m)
  colnames(design) <- make.names(colnames(design), unique = TRUE)
  expect_true("conditionbenign" %in% colnames(design))
  expect_false("conditioncontrol" %in% colnames(design))
})
