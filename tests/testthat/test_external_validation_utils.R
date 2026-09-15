test_that("standardize_gene_matrix collapses duplicated symbols", {
  x <- data.frame(gene = c("A", "A", "B"), s1 = c(1, 3, 2), s2 = c(2, 4, 3))
  mat <- standardize_gene_matrix(x)
  expect_true(all(c("A", "B") %in% rownames(mat)))
  expect_equal(nrow(mat), 2)
})

test_that("signed module scoring returns one value per sample", {
  x <- matrix(c(1,2,3,4, 4,3,2,1, 1,1,1,1), nrow = 3, byrow = TRUE)
  rownames(x) <- c("UP1", "DN1", "OTHER")
  colnames(x) <- paste0("s", 1:4)
  sc <- signed_module_score(x, up_genes = "UP1", down_genes = "DN1", min_genes = 1)
  expect_equal(length(sc), 4)
  expect_true(all(names(sc) == colnames(x)))
})
