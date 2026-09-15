# GEO acquisition -------------------------------------------------------------

select_eset <- function(gse_obj, platform_preference = NULL) {
  if (inherits(gse_obj, "ExpressionSet")) return(gse_obj)
  stopifnot(is.list(gse_obj), length(gse_obj) >= 1)
  if (!is.null(platform_preference)) {
    platforms <- vapply(gse_obj, function(es) Biobase::annotation(es) %||% "", character(1))
    idx <- which(platforms == platform_preference)
    if (length(idx)) return(gse_obj[[idx[[1]]]])
  }
  # Choose the expression set with the largest sample count.
  ns <- vapply(gse_obj, function(es) ncol(Biobase::exprs(es)), integer(1))
  gse_obj[[which.max(ns)]]
}

download_geo_one <- function(gse_id, ds_cfg, cfg) {
  cache <- file.path(cfg$project$data_dir, "raw", paste0(gse_id, "_eset.rds"))
  if (file.exists(cache)) return(readRDS(cache))
  message("Downloading ", gse_id, " from GEO...")
  g <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE, getGPL = TRUE, AnnotGPL = TRUE)
  es <- select_eset(g, ds_cfg$platform_preference %||% NULL)
  saveRDS(es, cache)
  es
}

download_all_geo <- function(cfg) {
  ids <- names(cfg$datasets)
  out <- stats::setNames(vector("list", length(ids)), ids)
  for (id in ids) out[[id]] <- download_geo_one(id, cfg$datasets[[id]], cfg)
  out
}
