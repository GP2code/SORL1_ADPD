#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

# ------------------------------------------------------------------
# Purpose:
#   Summarize SORL1 rare/damaging segregation per family and per variant
#   from PLINK2 --export A .raw matrices, joined to the informative PED.
#
# Output location:
#   By default, reads/writes next to this script (script-local).
#
# Env overrides:
#   OUT_DIR
#   PED_PATH
#   RAW_BASE_DIR
#   PVAR_BASE_DIR
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

ped_path <- Sys.getenv("PED_PATH", unset = file.path(out_dir, "r10_master_ped_informative_families.tsv"))
raw_base_dir <- Sys.getenv("RAW_BASE_DIR", unset = file.path(out_dir, "plink_chr11_SORL1_rare_damaging"))
pvar_base_dir <- Sys.getenv("PVAR_BASE_DIR", unset = raw_base_dir)

out_by_family  <- file.path(out_dir, "r10_SORL1_segregation_by_family.tsv")
out_by_variant <- file.path(out_dir, "r10_SORL1_segregation_by_variant.tsv")
out_qc         <- file.path(out_dir, "qc_SORL1_segregation_summary.tsv")

if (!file.exists(ped_path)) stop("PED_PATH not found: ", ped_path)
if (!dir.exists(raw_base_dir)) stop("RAW_BASE_DIR not found: ", raw_base_dir)

cat("OUT_DIR      :", out_dir, "\n")
cat("PED_PATH     :", ped_path, "\n")
cat("RAW_BASE_DIR :", raw_base_dir, "\n")
cat("PVAR_BASE_DIR:", pvar_base_dir, "\n\n")

# ------------------------------------------------------------------
# 1) Read informative master ped
# ------------------------------------------------------------------

ped <- read_tsv(ped_path, col_types = cols())

needed <- c("family_id", "id", "phenotype", "is_informative")
missing_cols <- setdiff(needed, names(ped))
if (length(missing_cols) > 0) {
  stop("Required columns missing in informative PED: ", paste(missing_cols, collapse = ", "))
}

# Use ancestry if present; otherwise fall back to wgs_label
if (!("ancestry" %in% names(ped)) && !("wgs_label" %in% names(ped))) {
  stop("Need ancestry labeling in PED: expected 'ancestry' and/or 'wgs_label'.")
}

ped_inf <- ped %>%
  filter(is_informative %in% TRUE) %>%
  mutate(
    phenotype_i = suppressWarnings(as.integer(phenotype)),
    ancestry_clean = if ("ancestry" %in% names(ped)) na_if(trimws(as.character(ancestry)), "") else NA_character_,
    wgs_label_clean = if ("wgs_label" %in% names(ped)) na_if(trimws(as.character(wgs_label)), "") else NA_character_,
    ancestry_final = coalesce(ancestry_clean, wgs_label_clean)
  ) %>%
  select(family_id, id, phenotype = phenotype_i, ancestry_final)

cat("Informative individuals in ped:", nrow(ped_inf), "\n")
cat("Missing ancestry among informative ped:", sum(is.na(ped_inf$ancestry_final) | ped_inf$ancestry_final == ""), "\n\n")

# ------------------------------------------------------------------
# 2) Auto-detect ancestries from available .raw files
# ------------------------------------------------------------------

raw_files <- list.files(raw_base_dir, pattern = "\\.raw$", recursive = TRUE, full.names = TRUE)
# Expect paths like: <RAW_BASE_DIR>/<ANC>/chr11_<ANC>_r10_SORL1_rare_damaging.raw
raw_tbl <- tibble(raw_path = raw_files) %>%
  mutate(
    ancestry = basename(dirname(raw_path)),
    fname = basename(raw_path)
  ) %>%
  filter(str_detect(fname, "_r10_SORL1_rare_damaging\\.raw$"))

ancestries <- sort(unique(raw_tbl$ancestry))
if (length(ancestries) == 0) {
  stop("No *.raw files found under RAW_BASE_DIR matching *_r10_SORL1_rare_damaging.raw")
}

cat("Ancestries detected from .raw files:", paste(ancestries, collapse = ", "), "\n\n")

# ------------------------------------------------------------------
# 3) Helper: robust read of PLINK .raw (tab or space)
# ------------------------------------------------------------------

read_plink_raw <- function(path) {
  # read_table2 handles both tabs and variable spaces well
  read_table2(path, col_types = cols())
}

# Convert dosage safely:
# - parse numeric
# - set values outside [0,2] (e.g., -9) to NA
to_dosage <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y[is.na(y)] <- NA_real_
  y[y < 0 | y > 2] <- NA_real_
  y
}

# ------------------------------------------------------------------
# 4) Process one ancestry
# ------------------------------------------------------------------

process_ancestry <- function(anc) {
  cat("==== Processing ancestry:", anc, "====\n")

  raw_path <- file.path(raw_base_dir, anc, paste0("chr11_", anc, "_r10_SORL1_rare_damaging.raw"))
  pvar_path <- file.path(pvar_base_dir, anc, paste0("chr11_", anc, "_r10_SORL1_rare_damaging.pvar"))

  if (!file.exists(raw_path)) {
    cat("  ! No .raw file at", raw_path, "– skipping.\n\n")
    return(list(by_family = tibble(), by_variant = tibble(), qc = tibble()))
  }
  if (!file.exists(pvar_path)) {
    cat("  ! No .pvar file at", pvar_path, "– skipping (needed for REF/ALT mapping).\n\n")
    return(list(by_family = tibble(), by_variant = tibble(), qc = tibble()))
  }

  raw <- read_plink_raw(raw_path)

  # Required ID columns in .raw
  req_raw_cols <- c("FID", "IID")
  if (!all(req_raw_cols %in% names(raw))) {
    cat("  ! .raw file missing FID/IID columns – skipping.\n\n")
    return(list(by_family = tibble(), by_variant = tibble(), qc = tibble()))
  }

  # Identify genotype columns
  fixed_cols <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
  geno_cols <- setdiff(names(raw), fixed_cols)

  if (length(geno_cols) == 0) {
    cat("  ! No genotype columns in .raw – skipping.\n\n")
    return(list(by_family = tibble(), by_variant = tibble(), qc = tibble()))
  }

  # Join with ped for this ancestry
  ped_anc <- ped_inf %>%
    filter(ancestry_final == anc) %>%
    select(family_id, id, phenotype)

  dat <- raw %>%
    rename(family_id = FID, id = IID) %>%
    inner_join(ped_anc, by = c("family_id", "id"))

  if (nrow(dat) == 0) {
    cat("  ! No overlap between .raw and ped for", anc, "– skipping.\n\n")
    return(list(by_family = tibble(), by_variant = tibble(), qc = tibble()))
  }

  cat("  Samples in .raw after join:", nrow(dat), "\n")
  cat("  Variants in .raw:", length(geno_cols), "\n")

  # Read pvar for REF/ALT mapping
  pvar <- read_tsv(pvar_path, col_types = cols(.default = col_character()))
  if (!all(c("#CHROM", "POS", "ID", "REF", "ALT") %in% names(pvar))) {
    cat("  ! pvar missing required columns (#CHROM, POS, ID, REF, ALT) – skipping.\n\n")
    return(list(by_family = tibble(), by_variant = tibble(), qc = tibble()))
  }

  id_map <- pvar %>%
    transmute(
      variant_id = ID,
      REF = REF,
      ALT = ALT
    ) %>%
    distinct(variant_id, .keep_all = TRUE)

  # Long format for summarization
  long <- dat %>%
    select(family_id, id, phenotype, all_of(geno_cols)) %>%
    pivot_longer(
      cols = all_of(geno_cols),
      names_to = "variant_col",
      values_to = "counted_dosage_raw"
    ) %>%
    mutate(
      # PLINK .raw columns are typically like "<variant_id>_<counted_allele>"
      variant_id = sub("_[^_]+$", "", variant_col),
      counted_allele = sub("^.*_", "", variant_col),
      counted_dosage = to_dosage(counted_dosage_raw)
    ) %>%
    select(-counted_dosage_raw) %>%
    left_join(id_map, by = "variant_id")

  # Compute ALT allele dosage robustly:
  # - if counted allele == ALT: alt_dosage = counted
  # - if counted allele == REF: alt_dosage = 2 - counted
  # - else: NA (allele mismatch)
  long <- long %>%
    mutate(
      alt_dosage = case_when(
        !is.na(ALT) & counted_allele == ALT ~ counted_dosage,
        !is.na(REF) & counted_allele == REF ~ if_else(is.na(counted_dosage), NA_real_, 2 - counted_dosage),
        TRUE ~ NA_real_
      ),
      is_nonmissing = !is.na(alt_dosage),
      is_carrier = is_nonmissing & alt_dosage > 0,
      is_PD   = phenotype == 2L,
      is_ctrl = phenotype == 1L,
      is_other = !(phenotype %in% c(1L, 2L))
    )

  # QC metrics
  n_missing_refalt <- sum(is.na(long$REF) | is.na(long$ALT))
  n_allele_mismatch <- sum(!(is.na(long$REF) | is.na(long$ALT)) & !(long$counted_allele %in% c(long$REF, long$ALT)))
  qc <- tibble(
    ancestry = anc,
    n_samples_joined = n_distinct(long$id),
    n_geno_columns = length(geno_cols),
    n_long_rows = nrow(long),
    n_rows_missing_REFALT = as.integer(n_missing_refalt),
    n_rows_allele_mismatch = as.integer(n_allele_mismatch),
    n_rows_missing_alt_dosage = as.integer(sum(!long$is_nonmissing))
  )

  if (n_allele_mismatch > 0) {
    cat("  NOTE: allele mismatches detected in .raw column suffix vs pvar REF/ALT. Rows affected:", n_allele_mismatch, "\n")
  }

  # Per-family, per-variant segregation
  by_family <- long %>%
    group_by(family_id, variant_col, variant_id) %>%
    summarise(
      ancestry          = anc,
      n_nonmissing      = sum(is_nonmissing),
      n_PD              = sum(is_nonmissing & is_PD),
      n_ctrl            = sum(is_nonmissing & is_ctrl),
      n_other           = sum(is_nonmissing & is_other),
      n_carrier_PD      = sum(is_carrier & is_PD),
      n_carrier_ctrl    = sum(is_carrier & is_ctrl),
      n_carrier_other   = sum(is_carrier & is_other),
      n_noncarrier_PD   = sum(is_nonmissing & !is_carrier & is_PD),
      n_noncarrier_ctrl = sum(is_nonmissing & !is_carrier & is_ctrl),
      n_noncarrier_other= sum(is_nonmissing & !is_carrier & is_other),
      .groups = "drop"
    ) %>%
    filter(n_nonmissing > 0) %>%
    arrange(family_id, variant_id)

  # Variant-level summary: collapse across families
  by_variant <- by_family %>%
    group_by(ancestry, variant_col, variant_id) %>%
    summarise(
      n_families        = n_distinct(family_id),
      n_nonmissing      = sum(n_nonmissing),
      n_PD              = sum(n_PD),
      n_ctrl            = sum(n_ctrl),
      n_other           = sum(n_other),
      n_carrier_PD      = sum(n_carrier_PD),
      n_carrier_ctrl    = sum(n_carrier_ctrl),
      n_carrier_other   = sum(n_carrier_other),
      n_noncarrier_PD   = sum(n_noncarrier_PD),
      n_noncarrier_ctrl = sum(n_noncarrier_ctrl),
      n_noncarrier_other= sum(n_noncarrier_other),
      .groups = "drop"
    ) %>%
    arrange(variant_id)

  cat("  Family-level rows :", nrow(by_family), "\n")
  cat("  Variant-level rows:", nrow(by_variant), "\n\n")

  list(by_family = by_family, by_variant = by_variant, qc = qc)
}

# ------------------------------------------------------------------
# 5) Run and write outputs
# ------------------------------------------------------------------

res_list <- lapply(ancestries, process_ancestry)

by_family_all  <- bind_rows(map(res_list, "by_family"))
by_variant_all <- bind_rows(map(res_list, "by_variant"))
qc_all         <- bind_rows(map(res_list, "qc")) %>% arrange(ancestry)

cat("Combined family-level rows :", nrow(by_family_all), "\n")
cat("Combined variant-level rows:", nrow(by_variant_all), "\n\n")

write_tsv(by_family_all,  out_by_family)
write_tsv(by_variant_all, out_by_variant)
write_tsv(qc_all,         out_qc)

cat("Wrote per-family segregation summary to:\n  ", out_by_family, "\n", sep = "")
cat("Wrote per-variant segregation summary to:\n  ", out_by_variant, "\n", sep = "")
cat("Wrote QC summary to:\n  ", out_qc, "\n", sep = "")

