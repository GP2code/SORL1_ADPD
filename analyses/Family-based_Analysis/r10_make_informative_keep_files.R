#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

# ------------------------------------------------------------------
# Purpose:
#   From the informative master PED, produce ancestry-specific PLINK keep
#   files (FID/IID) for informative WGS individuals.
#
# Output location:
#   By default, writes outputs next to this script (script-local).
#   Override with OUT_DIR if needed.
#
# Env overrides:
#   OUT_DIR
#   PED_INFORM_PATH
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

ped_inform_path <- Sys.getenv(
  "PED_INFORM_PATH",
  unset = file.path(out_dir, "r10_master_ped_informative_families.tsv")
)

if (!file.exists(ped_inform_path)) {
  stop(
    "Informative PED not found: ", ped_inform_path, "\n",
    "Run r10_mark_informative_families.R in the same folder first, or set PED_INFORM_PATH."
  )
}

qc_counts_path   <- file.path(out_dir, "qc_informative_wgs_counts_by_ancestry.tsv")
qc_missing_anc   <- file.path(out_dir, "qc_informative_wgs_missing_ancestry.tsv")
out_keep_all     <- file.path(out_dir, "ALL_r10_informative_keep.tsv")

cat("Output directory (script-local by default):\n  ", out_dir, "\n\n", sep = "")

# ------------------------------------------------------------------
# 1) Read informative master ped
# ------------------------------------------------------------------

ped <- read_tsv(ped_inform_path, col_types = cols())

cat("Total individuals in r10 master ped with informative flags:", nrow(ped), "\n")

required_cols <- c("family_id", "id", "has_wgs", "is_informative", "phenotype")
missing_cols <- setdiff(required_cols, names(ped))
if (length(missing_cols) > 0) {
  stop("Required columns missing from r10_master_ped_informative_families.tsv: ",
       paste(missing_cols, collapse = ", "))
}

# We prefer 'ancestry' if present; otherwise we will use master-key wgs_label if present.
if (!("ancestry" %in% names(ped)) && !("wgs_label" %in% names(ped))) {
  stop("Need ancestry labeling: expected column 'ancestry' and/or 'wgs_label' in informative PED.")
}

# ------------------------------------------------------------------
# 2) Restrict to informative WGS individuals (and build ancestry_final)
# ------------------------------------------------------------------

ped2 <- ped %>%
  mutate(
    # Robust phenotype coding
    phenotype_i = suppressWarnings(as.integer(phenotype)),

    # Robust ancestry label: use ancestry if present, else fall back to wgs_label
    ancestry_clean = if ("ancestry" %in% names(ped)) na_if(trimws(as.character(ancestry)), "") else NA_character_,
    wgs_label_clean = if ("wgs_label" %in% names(ped)) na_if(trimws(as.character(wgs_label)), "") else NA_character_,
    ancestry_final = dplyr::coalesce(ancestry_clean, wgs_label_clean)
  ) %>%
  select(-ancestry_clean, -wgs_label_clean)

inf_wgs <- ped2 %>%
  filter(is_informative %in% TRUE, has_wgs %in% TRUE)

cat("Informative WGS individuals (all):", nrow(inf_wgs), "\n")

# QC: missing ancestry among informative WGS
missing_anc_df <- inf_wgs %>%
  filter(is.na(ancestry_final) | ancestry_final == "") %>%
  transmute(
    family_id,
    id,
    phenotype = phenotype_i,
    ancestry = if ("ancestry" %in% names(ped2)) as.character(ancestry) else NA_character_,
    wgs_label = if ("wgs_label" %in% names(ped2)) as.character(wgs_label) else NA_character_
  )

write_tsv(missing_anc_df, qc_missing_anc)
cat("Informative WGS individuals with missing ancestry_final:", nrow(missing_anc_df), "\n")
cat("Wrote missing-ancestry QC to:\n  ", qc_missing_anc, "\n\n", sep = "")

# For keep-file generation, drop missing ancestry_final (but we keep QC above)
inf_wgs_anc <- inf_wgs %>%
  filter(!is.na(ancestry_final), ancestry_final != "")

cat("Informative WGS individuals (with ancestry_final):", nrow(inf_wgs_anc), "\n\n")

# ------------------------------------------------------------------
# 3) Summarize by ancestry (PD vs controls) and write QC table
# ------------------------------------------------------------------

summ_by_anc <- inf_wgs_anc %>%
  mutate(
    status = case_when(
      phenotype_i == 2L ~ "PD",
      phenotype_i == 1L ~ "Control",
      TRUE              ~ "Other"
    )
  ) %>%
  count(ancestry = ancestry_final, status, name = "n") %>%
  tidyr::pivot_wider(
    names_from  = status,
    values_from = n,
    values_fill = 0
  ) %>%
  arrange(ancestry)

write_tsv(summ_by_anc, qc_counts_path)

cat("Informative WGS counts by ancestry:\n")
print(summ_by_anc)
cat("\nWrote ancestry count QC to:\n  ", qc_counts_path, "\n\n", sep = "")

# ------------------------------------------------------------------
# 4) Write per-ancestry PLINK keep files (FID, IID) + combined keep
# ------------------------------------------------------------------

# Combined keep (all informative WGS with ancestry_final present)
keep_all <- inf_wgs_anc %>%
  transmute(FID = family_id, IID = id) %>%
  distinct()

write_tsv(keep_all, out_keep_all, col_names = FALSE)
cat("Wrote combined keep file (ALL ancestries) with", nrow(keep_all), "samples to:\n  ",
    out_keep_all, "\n\n", sep = "")

# Per-ancestry keep files
ancestries <- sort(unique(inf_wgs_anc$ancestry_final))

for (anc in ancestries) {
  sub <- inf_wgs_anc %>%
    filter(ancestry_final == anc) %>%
    transmute(FID = family_id, IID = id) %>%
    distinct()

  out_keep <- file.path(out_dir, paste0(anc, "_r10_informative_keep.tsv"))
  write_tsv(sub, out_keep, col_names = FALSE)

  cat("Wrote keep file for", anc, "with", nrow(sub), "samples to:\n  ", out_keep, "\n\n")
}

