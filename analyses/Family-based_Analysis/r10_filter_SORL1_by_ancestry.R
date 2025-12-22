#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# ------------------------------------------------------------------
# Purpose:
#   Filter ANNOVAR multianno tables to define "rare damaging" SORL1 variants
#   per ancestry, and output:
#     - filtered variant table (TSV)
#     - variant ID list for PLINK2 extraction
#
# Output location:
#   By default, reads/writes next to this script (script-local).
#
# Env overrides:
#   OUT_DIR
#   ANNOVAR_DIR
#   PFILE_BASE_DIR
#   AF_CUTOFF (default 0.01)
#   CADD_CUTOFF (default 20)
# ------------------------------------------------------------------

Sys.setenv(TZ = "UTC")
options(readr.show_col_types = FALSE)

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 0L) return(getwd())
  script_path <- sub("^--file=", "", file_arg[1])
  dirname(normalizePath(script_path))
}

script_dir <- get_script_dir()
out_dir <- Sys.getenv("OUT_DIR", unset = script_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ANNOVAR_DIR <- Sys.getenv("ANNOVAR_DIR", unset = file.path(out_dir, "annovar_r10"))
PFILE_BASE_DIR <- Sys.getenv("PFILE_BASE_DIR", unset = file.path(out_dir, "plink_chr11_informative"))

AF_CUTOFF <- suppressWarnings(as.numeric(Sys.getenv("AF_CUTOFF", unset = "0.01")))
if (is.na(AF_CUTOFF)) AF_CUTOFF <- 0.01

CADD_CUTOFF <- suppressWarnings(as.numeric(Sys.getenv("CADD_CUTOFF", unset = "20")))
if (is.na(CADD_CUTOFF)) CADD_CUTOFF <- 20

if (!dir.exists(ANNOVAR_DIR)) stop("ANNOVAR_DIR not found: ", ANNOVAR_DIR)

qc_summary_path <- file.path(out_dir, "qc_SORL1_filter_summary.tsv")
qc_rows <- list()

cat("OUT_DIR     :", out_dir, "\n")
cat("ANNOVAR_DIR :", ANNOVAR_DIR, "\n")
cat("PFILE_BASE_DIR (for ID mapping):", PFILE_BASE_DIR, "\n")
cat("AF_CUTOFF   :", AF_CUTOFF, "\n")
cat("CADD_CUTOFF :", CADD_CUTOFF, "\n\n")

# ------------------------------------------------------------------
# 1) Map GP2 ancestry -> gnomAD subpop tags (lowercase)
# ------------------------------------------------------------------

gnomad_tags_by_anc <- list(
  AFR = c("afr"),
  AJ  = c("asj"),
  AMR = c("amr", "lat"),
  CAH = character(0),    # use global + popmax + grpmax only
  CAS = character(0),    # use global + popmax + grpmax only
  EAS = c("eas"),
  EUR = c("nfe", "eur"),
  FIN = c("fin"),
  MDE = c("mid"),
  SAS = c("sas")
)

# ------------------------------------------------------------------
# 2) Helper: pick AF columns for this ancestry (global + popmax + grpmax + subpop)
# ------------------------------------------------------------------

pick_af_cols <- function(anc, colnames_vec) {
  af_cols_all <- grep("^gnomad41_(exome|genome)_AF", colnames_vec, value = TRUE)

  global_af <- af_cols_all[grepl("_AF$", af_cols_all)]
  popmax_af <- af_cols_all[grepl("_AF_popmax$", af_cols_all)]
  grpmax_af <- af_cols_all[grepl("_AF_grpmax$", af_cols_all)]

  subpop_af <- af_cols_all[
    grepl("_AF_", af_cols_all) &
      !grepl("_AF_popmax$", af_cols_all) &
      !grepl("_AF_grpmax$", af_cols_all)
  ]

  tags <- gnomad_tags_by_anc[[anc]]
  if (is.null(tags)) tags <- character(0)

  subpop_sel <- character(0)
  if (length(tags) > 0 && length(subpop_af) > 0) {
    subpop_sel <- subpop_af[
      vapply(subpop_af, function(col) {
        suf <- tolower(sub(".*_AF_", "", col))
        any(suf %in% tags)
      }, logical(1))
    ]
  }

  unique(c(global_af, popmax_af, grpmax_af, subpop_sel))
}

# ------------------------------------------------------------------
# 3) LOF definition (no stoploss)
# ------------------------------------------------------------------

is_lof_variant_vec <- function(exonic_func, func_ref) {
  # Vectorized LOF logic, tolerant of compound annotations
  lof_exonic <- c(
    "stopgain",
    "frameshift insertion",
    "frameshift deletion",
    "frameshift block substitution"
  )
  exonic_is_lof <- vapply(exonic_func, function(x) {
    if (is.na(x)) return(FALSE)
    any(str_detect(x, fixed(lof_exonic)))
  }, logical(1))

  func_is_splicing <- !is.na(func_ref) & str_detect(func_ref, "splicing")

  exonic_is_lof | func_is_splicing
}

# Robust numeric conversion (treat "." and "" as NA)
to_num <- function(x) suppressWarnings(as.numeric(na_if(na_if(x, "."), "")))

# Robust SORL1 token match in Gene.refGene
is_sorl1_gene <- function(gene_field) {
  gf <- coalesce(gene_field, "")
  # Split tokens on common separators; match exact "SORL1"
  tokens <- str_split(gf, pattern = "[,;]", simplify = FALSE)
  vapply(tokens, function(toks) {
    toks2 <- trimws(toks)
    any(toks2 == "SORL1")
  }, logical(1))
}

# ------------------------------------------------------------------
# 4) Loop over ancestries present in annovar_r10/
# ------------------------------------------------------------------

anc_dirs <- list.dirs(ANNOVAR_DIR, full.names = FALSE, recursive = FALSE)
anc_dirs <- anc_dirs[anc_dirs != ""]

cat("Ancestries found in ANNOVAR_DIR:\n")
print(anc_dirs)
cat("\n")

for (anc in anc_dirs) {
  cat("==== Processing ancestry:", anc, "====\n")

  anc_dir <- file.path(ANNOVAR_DIR, anc)

  # Prefer the expected file name if present; otherwise take the first multianno.
  expected_prefix <- file.path(anc_dir, paste0("chr11_", anc, "_r10_informative_SORL1.annotated.hg38_multianno.txt"))
  multianno_files <- list.files(anc_dir, pattern = "hg38_multianno\\.txt$", full.names = TRUE)

  if (length(multianno_files) == 0) {
    cat("  ! No multianno files in", anc_dir, "– skipping.\n\n")
    next
  }

  multianno_path <- if (file.exists(expected_prefix)) expected_prefix else multianno_files[1]
  if (!file.exists(multianno_path)) {
    cat("  ! Selected multianno file does not exist:", multianno_path, "\n\n")
    next
  }

  if (length(multianno_files) > 1 && !file.exists(expected_prefix)) {
    cat("  NOTE: Multiple multianno files found; using first:", multianno_path, "\n")
  } else {
    cat("  Using multianno file:", multianno_path, "\n")
  }

  ann <- read_tsv(multianno_path, col_types = cols(.default = col_character()))
  colnames_ann <- names(ann)

  required_cols <- c("Chr", "Start", "End", "Ref", "Alt",
                     "Gene.refGene", "Func.refGene", "ExonicFunc.refGene")
  missing_cols <- setdiff(required_cols, colnames_ann)
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", multianno_path, ": ", paste(missing_cols, collapse = ", "))
  }

  cadd_col <- "CADD_phred"
  if (!cadd_col %in% colnames_ann) {
    stop("Could not find CADD_phred column in ", multianno_path,
         ". If your dbNSFP CADD column name differs, update cadd_col.")
  }

  # Restrict to SORL1 (robust token match)
  keep_sorl1 <- is_sorl1_gene(ann[["Gene.refGene"]])
  ann_sorl1 <- ann[keep_sorl1, , drop = FALSE]

  n_total_sorl1 <- nrow(ann_sorl1)
  cat("  Total SORL1 rows:", n_total_sorl1, "\n")

  if (n_total_sorl1 == 0) {
    qc_rows[[length(qc_rows) + 1L]] <- tibble(
      ancestry = anc,
      n_sorl1_rows = 0L,
      n_filtered = 0L,
      n_mapped_to_pvar_id = 0L,
      n_unmapped = 0L,
      af_cols_used = "",
      multianno = basename(multianno_path)
    )
    cat("  No SORL1 variants for", anc, "– skipping.\n\n")
    next
  }

  # AF columns for this ancestry
  af_cols <- pick_af_cols(anc, colnames_ann)
  cat("  AF columns used for", anc, ":", if (length(af_cols) == 0) "<none>" else paste(af_cols, collapse = ", "), "\n")

  # Compute numeric columns
  CADD_phred_num <- to_num(ann_sorl1[[cadd_col]])

  if (length(af_cols) > 0) {
    af_mat <- sapply(af_cols, function(col) to_num(ann_sorl1[[col]]))
    if (!is.matrix(af_mat)) af_mat <- matrix(af_mat, ncol = length(af_cols))
    max_AF_vec <- apply(af_mat, 1, function(r) if (all(is.na(r))) NA_real_ else max(r, na.rm = TRUE))
  } else {
    max_AF_vec <- rep(NA_real_, nrow(ann_sorl1))
  }

  # Functional class flags
  is_LOF <- is_lof_variant_vec(ann_sorl1[["ExonicFunc.refGene"]], ann_sorl1[["Func.refGene"]])
  is_mis <- !is.na(ann_sorl1[["ExonicFunc.refGene"]]) & str_detect(ann_sorl1[["ExonicFunc.refGene"]], "nonsynonymous SNV")

  ann_sorl1_num <- ann_sorl1 %>%
    mutate(
      CADD_phred_num = CADD_phred_num,
      max_AF = max_AF_vec,
      is_LOF = is_LOF,
      is_mis = is_mis
    )

  # Rare + damaging filter
  filtered <- ann_sorl1_num %>%
    filter(
      (is_LOF | (is_mis & !is.na(CADD_phred_num) & CADD_phred_num >= CADD_CUTOFF)),
      is.na(max_AF) | max_AF <= AF_CUTOFF
    )

  n_filtered <- nrow(filtered)
  cat("  Rare/damaging SORL1 variants passing filters:", n_filtered, "\n")

  if (n_filtered == 0) {
    qc_rows[[length(qc_rows) + 1L]] <- tibble(
      ancestry = anc,
      n_sorl1_rows = as.integer(n_total_sorl1),
      n_filtered = 0L,
      n_mapped_to_pvar_id = 0L,
      n_unmapped = 0L,
      af_cols_used = paste(af_cols, collapse = ","),
      multianno = basename(multianno_path)
    )
    cat("\n")
    next
  }

  # Build clean key for mapping: CHR:POS:REF:ALT with chr stripped
  chr_clean <- str_replace(tolower(filtered$Chr), "^chr", "")
  pos_clean <- filtered$Start
  ref_clean <- filtered$Ref
  alt_clean <- filtered$Alt
  key_clean <- paste0(chr_clean, ":", pos_clean, ":", ref_clean, ":", alt_clean)

  # Try to map to actual PLINK2 variant IDs from the extracted pvar (best practice)
  pvar_path <- file.path(PFILE_BASE_DIR, anc, paste0("chr11_", anc, "_r10_informative_SORL1.pvar"))
  id_map <- NULL
  if (file.exists(pvar_path)) {
    pvar <- read_tsv(pvar_path, col_types = cols(.default = col_character()))
    if (all(c("#CHROM", "POS", "ID", "REF", "ALT") %in% names(pvar))) {
      p_chr <- str_replace(tolower(pvar[["#CHROM"]]), "^chr", "")
      p_key <- paste0(p_chr, ":", pvar[["POS"]], ":", pvar[["REF"]], ":", pvar[["ALT"]])
      id_map <- tibble(key_clean = p_key, plink_id = pvar[["ID"]]) %>%
        distinct(key_clean, .keep_all = TRUE)
    } else {
      cat("  NOTE: pvar did not have expected columns (#CHROM, POS, ID, REF, ALT); skipping ID mapping.\n")
    }
  } else {
    cat("  NOTE: pvar not found for ID mapping (", pvar_path, "); will output computed IDs.\n", sep = "")
  }

  filtered2 <- filtered %>%
    mutate(key_clean = key_clean)

  if (!is.null(id_map)) {
    filtered2 <- filtered2 %>%
      left_join(id_map, by = "key_clean")
  } else {
    filtered2 <- filtered2 %>%
      mutate(plink_id = NA_character_)
  }

  # Final variant_id: prefer PLINK ID, else fall back to computed "chr<Chr>:<Start>:<Ref>:<Alt>"
  filtered2 <- filtered2 %>%
    mutate(
      variant_id = if_else(
        !is.na(plink_id) & plink_id != "",
        plink_id,
        paste0("chr", str_replace(tolower(Chr), "^chr", ""), ":", Start, ":", Ref, ":", Alt)
      )
    )

  n_mapped <- sum(!is.na(filtered2$plink_id) & filtered2$plink_id != "")
  n_unmapped <- n_filtered - n_mapped
  cat("  ID mapping: mapped", n_mapped, " / ", n_filtered, " variants to pvar IDs.\n\n", sep = "")

  # Write outputs
  out_tsv <- file.path(anc_dir, paste0("chr11_", anc, "_r10_SORL1_rare_damaging_variants.tsv"))
  out_ids <- file.path(anc_dir, paste0("chr11_", anc, "_r10_SORL1_rare_damaging_ids.txt"))

  write_tsv(filtered2, out_tsv)
  writeLines(filtered2$variant_id, out_ids)

  cat("  Wrote filtered variant table to:", out_tsv, "\n")
  cat("  Wrote variant ID list to     :", out_ids, "\n\n")

  qc_rows[[length(qc_rows) + 1L]] <- tibble(
    ancestry = anc,
    n_sorl1_rows = as.integer(n_total_sorl1),
    n_filtered = as.integer(n_filtered),
    n_mapped_to_pvar_id = as.integer(n_mapped),
    n_unmapped = as.integer(n_unmapped),
    af_cols_used = paste(af_cols, collapse = ","),
    multianno = basename(multianno_path)
  )
}

# Write QC summary across ancestries
qc_df <- if (length(qc_rows) == 0) {
  tibble(
    ancestry = character(),
    n_sorl1_rows = integer(),
    n_filtered = integer(),
    n_mapped_to_pvar_id = integer(),
    n_unmapped = integer(),
    af_cols_used = character(),
    multianno = character()
  )
} else {
  bind_rows(qc_rows) %>% arrange(ancestry)
}

write_tsv(qc_df, qc_summary_path)
cat("Wrote QC summary to:\n  ", qc_summary_path, "\n", sep = "")
cat("All ancestries processed for SORL1 rare/damaging filtering.\n")

