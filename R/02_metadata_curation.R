# Metadata curation -----------------------------------------------------------

empty_metadata_overrides <- function() {
  tibble::tibble(
    dataset = character(),
    sample = character(),
    condition = character(),
    patient_id = character(),
    technical_group = character(),
    include = logical(),
    notes = character()
  )
}

load_manual_overrides <- function(path = "config/manual_overrides/sample_metadata_overrides.csv") {
  # Keep a stable empty schema. A zero-row tibble without `dataset` caused
  # targets::tar_make() to fail in metadata_list when the manual override file
  # contained only the header/comment template.
  if (!file.exists(path)) return(empty_metadata_overrides())

  suppressWarnings({
    x <- readr::read_csv(path, comment = "#", show_col_types = FALSE)
  })
  x <- janitor::clean_names(x)

  required <- names(empty_metadata_overrides())
  for (cc in setdiff(required, names(x))) {
    x[[cc]] <- empty_metadata_overrides()[[cc]]
  }
  x <- x |> dplyr::select(dplyr::all_of(required))

  if (!nrow(x)) return(empty_metadata_overrides())

  x |>
    dplyr::mutate(
      dataset = as.character(.data$dataset),
      sample = as.character(.data$sample),
      condition = as.character(.data$condition),
      patient_id = as.character(.data$patient_id),
      technical_group = as.character(.data$technical_group),
      include = dplyr::case_when(
        is.na(.data$include) ~ NA,
        as.character(.data$include) %in% c("TRUE", "True", "true", "1", "yes", "YES") ~ TRUE,
        as.character(.data$include) %in% c("FALSE", "False", "false", "0", "no", "NO") ~ FALSE,
        TRUE ~ as.logical(.data$include)
      ),
      notes = as.character(.data$notes)
    )
}

infer_condition <- function(text, condition_regex) {
  text <- tolower(text)
  tumor <- grepl(condition_regex$tumor, text, ignore.case = TRUE, perl = TRUE)
  benign <- if (!is.null(condition_regex$benign)) grepl(condition_regex$benign, text, ignore.case = TRUE, perl = TRUE) else rep(FALSE, length(text))
  cp <- if (!is.null(condition_regex$cp)) grepl(condition_regex$cp, text, ignore.case = TRUE, perl = TRUE) else rep(FALSE, length(text))
  control <- grepl(condition_regex$control, text, ignore.case = TRUE, perl = TRUE)

  # Order matters. Many GEO records describe normal/adjacent samples as coming
  # from pancreatic-cancer patients, and terms such as "nontumor" contain the
  # substring "tumor". Therefore, explicit control/benign evidence must override
  # generic tumor/cancer evidence. This was the cause of the v1.1 one-level
  # contrast failure for GSE15471 and would also affect GSE28735.
  out <- rep(NA_character_, length(text))
  out[tumor] <- "tumor"
  out[benign] <- "benign"
  # Chronic pancreatitis is biologically distinct from benign neoplasm/adjacent benign tissue.
  # It should remain visible in metadata audits but should not be coerced into the
  # normal-benign-PDAC transition axis.
  out[cp] <- "cp"
  out[control] <- "control"
  out
}

metadata_biology_text <- function(pdata) {
  # GEO pData contains long protocol fields. Those fields often contain words
  # such as "hybridisation controls", "normalization", or "cancer patients"
  # that are not sample labels. For condition inference, use only columns likely
  # to carry biological/sample annotation, then dataset-specific rules below.
  nm <- names(pdata)
  keep <- grepl(
    paste(
      c(
        "^sample$", "title", "source", "characteristics", "tissue",
        "disease", "diagnosis", "status", "group", "sample_type",
        "sampletype", "phenotype", "description"
      ),
      collapse = "|"
    ),
    nm, ignore.case = TRUE
  )
  if (!any(keep)) keep <- rep(TRUE, length(nm))
  apply(as.data.frame(pdata[, keep, drop = FALSE]), 1, paste, collapse = " | ")
}

infer_condition_from_title <- function(gse_id, pdata, current_condition) {
  title_col <- first_present_col(pdata, c("title", "sample_title", "source_name_ch1", "characteristics_ch1"))
  title <- if (!is.na(title_col)) as.character(pdata[[title_col]]) else rep("", nrow(pdata))
  title_l <- tolower(title)
  bio_text <- tolower(metadata_biology_text(pdata))

  # GSE15471 uses compact titles such as N30162, N30162_rep, T30162, T30162_rep.
  if (identical(gse_id, "GSE15471")) {
    current_condition[grepl("^\\s*n[0-9]+(_rep)?\\s*$", title_l, perl = TRUE)] <- "control"
    current_condition[grepl("^\\s*t[0-9]+(_rep)?\\s*$", title_l, perl = TRUE)] <- "tumor"
  }

  # GSE28735 titles are explicit: "human pancreatic tumor tissue..." or
  # "human pancreatic nontumor tissue...". Match nontumor before tumor.
  if (identical(gse_id, "GSE28735")) {
    current_condition[grepl("non[- ]?tumou?r|nontumou?r|adjacent|tissue:\\s*n\\b", bio_text, perl = TRUE)] <- "control"
    current_condition[grepl("(^| )tumou?r tissue|pancreatic tumou?r tissue|tissue:\\s*t\\b", bio_text, perl = TRUE) &
                        !grepl("non[- ]?tumou?r|nontumou?r", bio_text, perl = TRUE)] <- "tumor"
  }

  # GSE62165 has protocol text containing "hybridisation controls" for every
  # sample. Use explicit sample fields instead: tissue/pdac/control status.
  if (identical(gse_id, "GSE62165")) {
    current_condition[grepl("tissue:\\s*non[- ]?tumou?ral pancreatic tissue|(^|\\|)\\s*control\\s*(\\||$)|control sample", bio_text, perl = TRUE)] <- "control"
    current_condition[grepl("tissue:\\s*pancreatic tumou?r|(^|\\|)\\s*pdac\\s*(\\||$)", bio_text, perl = TRUE)] <- "tumor"
  }

  # GSE16515 titles/source names are explicit: Pancreatic Sample X-Tumor/Normal.
  if (identical(gse_id, "GSE16515")) {
    # Use direct sample-title endings before any generic regex. Avoid raw GEO
    # protocol fields because they can contain words like "normalization".
    current_condition[grepl("[- ]normal\\s*$|normal tissue", title_l, perl = TRUE)] <- "control"
    current_condition[grepl("[- ]tumou?r\\s*$|tumou?r tissue", title_l, perl = TRUE)] <- "tumor"
  }

  # GSE71989 titles/status separate normal pancreatic tissue from PDAC tissue.
  if (identical(gse_id, "GSE71989")) {
    current_condition[grepl("subject status:\\s*normal|human normal pancreatic tissue|normal pancreatic tissue|^normal$", bio_text, perl = TRUE) |
                        grepl("\\bnormal\\b", title_l, perl = TRUE)] <- "control"
    current_condition[grepl("subject status:\\s*pdac|human pdac tissue|pdac tissue|pancreatic ductal adenocarcinoma|^pdac$", bio_text, perl = TRUE) |
                        grepl("\\bpdac\\b|ductal adenocarcinoma|tumou?r", title_l, perl = TRUE)] <- "tumor"
  }

  # GSE91035 is the transition cohort: normal -> benign -> PDAC.
  # Manual GEO audit confirmed that GSM2420007 and GSM2420010 are chronic
  # pancreatitis (CP), not benign. Keep CP as an explicit label so it is visible
  # in the audit, but exclude it from transition contrasts downstream.
  if (identical(gse_id, "GSE91035")) {
    current_condition[grepl("disease state:\\s*normal|human normal pancreatic tissue|normal pancreatic tissue", bio_text, perl = TRUE)] <- "control"
    current_condition[grepl("disease state:\\s*benign|human benign pancreatic tissue|benign pancreatic tissue", bio_text, perl = TRUE)] <- "benign"
    current_condition[grepl("disease state:\\s*chronic pancreatitis|chronic pancreatitis", bio_text, perl = TRUE)] <- "cp"
    current_condition[grepl("disease state:\\s*pdac|human pdac tissue|pdac tissue|pancreatic ductal adenocarcinoma", bio_text, perl = TRUE)] <- "tumor"
    current_condition[pdata$sample %in% c("GSM2420007", "GSM2420010")] <- "cp"
  }

  current_condition
}

infer_patient_id <- function(pdata, regex, gse_id = NULL) {
  n <- nrow(pdata)
  all_text <- apply(as.data.frame(pdata), 1, paste, collapse = " | ")
  pid <- extract_regex_last_group(all_text, regex)

  # Common GEO columns if available.
  candidate_cols <- grep("patient|case|subject|donor|individual|specimen|source_name|title", names(pdata), ignore.case = TRUE, value = TRUE)
  for (cc in candidate_cols) {
    tmp <- extract_regex_last_group(as.character(pdata[[cc]]), regex)
    pid[is.na(pid)] <- tmp[is.na(pid)]
  }

  # GSE15471 patient IDs are encoded in sample titles like N30162/T30162 and
  # technical replicate titles like N30162_rep/T30162_rep.
  if (identical(gse_id, "GSE15471")) {
    title_col <- first_present_col(pdata, c("title", "sample_title", "source_name_ch1"))
    if (!is.na(title_col)) {
      tmp <- extract_regex_last_group(as.character(pdata[[title_col]]), "^[NTnt]([0-9]+)(?:_rep)?$")
      pid[is.na(pid)] <- tmp[is.na(pid)]
    }
  }

  pid_num <- gsub("[^0-9]", "", as.character(pid))
  ifelse(
    is.na(pid) | !nzchar(pid_num),
    NA_character_,
    paste0("P", sprintf("%03d", as.integer(pid_num)))
  )
}

curate_metadata_one <- function(gse_id, es, ds_cfg, overrides = tibble::tibble()) {
  pdata <- Biobase::pData(es) |> as.data.frame() |> tibble::rownames_to_column("sample") |> janitor::clean_names()
  text <- apply(pdata, 1, paste, collapse = " | ")
  bio_text <- metadata_biology_text(pdata)
  title_value <- if ("title" %in% names(pdata)) as.character(pdata$title) else rep(NA_character_, nrow(pdata))

  meta <- pdata |>
    dplyr::transmute(
      dataset = gse_id,
      sample = .data$sample,
      title = title_value,
      raw_text = text,
      condition = infer_condition(bio_text, ds_cfg$condition_regex),
      patient_id = infer_patient_id(pdata, ds_cfg$patient_regex %||% "(patient|case|pt)[_ -]*0*([0-9]+)", gse_id = gse_id),
      include = TRUE
    )

  meta$condition <- infer_condition_from_title(gse_id, pdata, meta$condition)

  meta$technical_group <- NA_character_
  meta$notes <- NA_character_
  ov <- empty_metadata_overrides()
  if ("dataset" %in% names(overrides)) {
    ov <- overrides |> dplyr::filter(.data$dataset == gse_id)
  }
  if (nrow(ov)) {
    ov <- ov |> dplyr::select(dplyr::any_of(c("sample", "condition", "patient_id", "technical_group", "include", "notes")))
    for (i in seq_len(nrow(ov))) {
      idx <- match(ov$sample[[i]], meta$sample)
      if (is.na(idx)) next
      for (cc in intersect(names(ov), c("condition", "patient_id", "technical_group", "include", "notes"))) {
        val <- ov[[cc]][[i]]
        if (!is.na(val) && nzchar(as.character(val))) meta[[cc]][[idx]] <- val
      }
    }
    meta$include <- as.logical(meta$include)
  }

  # If patient IDs cannot be inferred for paired datasets, make a transparent warning.
  if (isTRUE(ds_cfg$paired_design) || identical(ds_cfg$paired_design, "mixed")) {
    if (sum(!is.na(meta$patient_id)) < 2) {
      warning(gse_id, ": patient_id could not be robustly inferred. Edit config/manual_overrides/sample_metadata_overrides.csv before final analysis.")
    }
  }

  meta <- meta |>
    dplyr::mutate(
      condition = dplyr::case_when(
        condition %in% c("pdac", "cancer") ~ "tumor",
        condition %in% c("normal", "adjacent", "control") ~ "control",
        condition %in% c("chronic pancreatitis", "pancreatitis", "cp") ~ "cp",
        TRUE ~ condition
      ),
      technical_group = dplyr::if_else(
        !is.na(.data$technical_group) & nzchar(.data$technical_group),
        .data$technical_group,
        paste(.data$dataset, dplyr::coalesce(.data$patient_id, .data$sample), .data$condition, sep = "__")
      )
    )

  valid_conditions <- c("tumor", "control", "benign", "cp")
  missing_condition <- meta |> dplyr::filter(is.na(.data$condition) | !.data$condition %in% valid_conditions)
  if (nrow(missing_condition)) {
    warning(gse_id, ": ", nrow(missing_condition), " samples lack a confident condition label; they will be excluded unless manually overridden.")
    meta$include[match(missing_condition$sample, meta$sample)] <- FALSE
  }

  meta
}

curate_all_metadata <- function(esets, cfg) {
  overrides <- load_manual_overrides()
  missing_cfg <- setdiff(names(esets), names(cfg$datasets))
  if (length(missing_cfg)) {
    stop("Missing dataset configuration for: ", paste(missing_cfg, collapse = ", "), call. = FALSE)
  }
  out <- purrr::imap(esets, ~ curate_metadata_one(.y, .x, cfg$datasets[[.y]], overrides))
  all <- dplyr::bind_rows(out)
  if (!"dataset" %in% names(all)) {
    stop("Metadata curation returned no `dataset` column. Check GEO download and metadata parsing.", call. = FALSE)
  }
  write_tsv(all, file.path(cfg$project$output_dir, "tables", "curated_metadata_all.tsv"))
  out
}

export_metadata_audit <- function(esets, cfg) {
  m <- curate_all_metadata(esets, cfg)
  invisible(m)
}
