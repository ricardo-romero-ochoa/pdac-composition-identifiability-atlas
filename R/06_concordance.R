# Cross-platform concordance --------------------------------------------------

pairwise_logfc_concordance <- function(de_list, cfg) {
  wide <- dplyr::bind_rows(de_list) |>
    dplyr::select(.data$dataset, .data$gene, .data$logFC) |>
    tidyr::pivot_wider(names_from = "dataset", values_from = "logFC")
  mat <- wide |> tibble::column_to_rownames("gene") |> as.matrix()
  pear <- stats::cor(mat, use = "pairwise.complete.obs", method = "pearson")
  spear <- stats::cor(mat, use = "pairwise.complete.obs", method = "spearman")
  write_tsv(as.data.frame(as.table(pear)) |> dplyr::rename(dataset1 = Var1, dataset2 = Var2, pearson = Freq), file.path(cfg$project$output_dir, "tables", "cross_dataset_logfc_pearson.tsv"))
  write_tsv(as.data.frame(as.table(spear)) |> dplyr::rename(dataset1 = Var1, dataset2 = Var2, spearman = Freq), file.path(cfg$project$output_dir, "tables", "cross_dataset_logfc_spearman.tsv"))
  list(wide = wide, pearson = pear, spearman = spear)
}

rank_overlap_one <- function(de_a, de_b, n = 500) {
  a_up <- de_a |> dplyr::arrange(dplyr::desc(.data$logFC)) |> dplyr::slice_head(n = n) |> dplyr::pull(.data$gene)
  b_up <- de_b |> dplyr::arrange(dplyr::desc(.data$logFC)) |> dplyr::slice_head(n = n) |> dplyr::pull(.data$gene)
  a_dn <- de_a |> dplyr::arrange(.data$logFC) |> dplyr::slice_head(n = n) |> dplyr::pull(.data$gene)
  b_dn <- de_b |> dplyr::arrange(.data$logFC) |> dplyr::slice_head(n = n) |> dplyr::pull(.data$gene)
  universe <- length(unique(c(de_a$gene, de_b$gene)))
  tibble::tibble(
    top_n = n,
    up_overlap = length(intersect(a_up, b_up)),
    up_jaccard = length(intersect(a_up, b_up)) / length(union(a_up, b_up)),
    up_p = stats::phyper(length(intersect(a_up, b_up)) - 1, n, universe - n, n, lower.tail = FALSE),
    down_overlap = length(intersect(a_dn, b_dn)),
    down_jaccard = length(intersect(a_dn, b_dn)) / length(union(a_dn, b_dn)),
    down_p = stats::phyper(length(intersect(a_dn, b_dn)) - 1, n, universe - n, n, lower.tail = FALSE)
  )
}

rank_overlap_concordance <- function(de_list, cfg) {
  ids <- names(de_list)
  ns <- unlist(cfg$analysis$top_n_rank_overlap %||% c(100, 250, 500, 1000))
  out <- purrr::map_dfr(utils::combn(ids, 2, simplify = FALSE), function(pair) {
    purrr::map_dfr(ns, function(n) rank_overlap_one(de_list[[pair[1]]], de_list[[pair[2]]], n = n)) |>
      dplyr::mutate(dataset1 = pair[1], dataset2 = pair[2], .before = 1)
  }) |>
    dplyr::mutate(up_fdr = bh(.data$up_p), down_fdr = bh(.data$down_p))
  write_tsv(out, file.path(cfg$project$output_dir, "tables", "cross_dataset_rank_overlap.tsv"))
  out
}

run_concordance <- function(de_list, cfg) {
  list(
    logfc = pairwise_logfc_concordance(de_list, cfg),
    rank_overlap = rank_overlap_concordance(de_list, cfg)
  )
}
