#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ------------------------------------------------------------------
# Purpose:
#   Annotate the r10 "core" pedigree with Release 10 master key fields
#   plus selected extended clinical age fields.
#
# Output location:
#   By default, writes output next to this script (script-local).
#   You can override by setting OUT_DIR environment variable.
#
# Env overrides:
#   CORE_PED_PATH, MASTER_KEY_PATH, EXT_CLIN_PATH, OUT_DIR
# ------------------------------------------------------------------

# Reduce noisy environment-related warnings/messages in container setups
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

# ------------------------------------------------------------------
# Paths (defaults)
# ------------------------------------------------------------------

core_ped_path <- Sys.getenv(
  "CORE_PED_PATH",
  unset = "/home/jupyter/workspace/ws_files/bernabe/analyses/r10_family_analysis/segregation_report_ped_r10.tsv"
)

master_key_path <- Sys.getenv(
  "MASTER_KEY_PATH",
  unset = "/home/jupyter/workspace/gp2_tier2_eu_release10/clinical_data/master_key_release10_final_vwb.csv"
)

ext_clin_path <- Sys.getenv(
  "EXT_CLIN_PATH",
  unset = "/home/jupyter/workspace/gp2_tier2_eu_release10/clinical_data/r10_extended_clinical_data_vwb.csv"
)

for (p in c(core_ped_path, master_key_path, ext_clin_path)) {
  if (!file.exists(p)) stop("Input file not found: ", p)
}

out_ped_annot <- file.path(out_dir, "r10_core_ped_with_master_ages.tsv")
qc_ext_dups   <- file.path(out_dir, "qc_extended_clinical_gp2id_duplicates.tsv")
qc_ext_conf   <- file.path(out_dir, "qc_extended_clinical_gp2id_conflicts.tsv")

cat("Output directory (script-local by default):\n  ", out_dir, "\n\n", sep = "")

# ------------------------------------------------------------------
# 1) Read r10 core family ped
# ------------------------------------------------------------------

core_ped <- read_tsv(
  core_ped_path,
  col_types = cols(
    family_id   = col_character(),
    id          = col_character(),
    parental_id = col_character(),
    maternal_id = col_character(),
    sex         = col_double(),
    phenotype   = col_double()
  )
)

dup_ids <- core_ped$id[duplicated(core_ped$id)]
if (length(dup_ids) > 0) {
  stop(
    "Duplicate 'id' values detected in core pedigree (would break joins). ",
    "Examples: ", paste(unique(dup_ids)[seq_len(min(10, length(unique(dup_ids))))], collapse = ", ")
  )
}

cat("Core r10 ped:\n")
cat("  Individuals:", nrow(core_ped), "\n")
cat("  Families   :", n_distinct(core_ped$family_id), "\n\n")

# ------------------------------------------------------------------
# 2) Read r10 master key and extended clinical data
# ------------------------------------------------------------------

mk <- suppressMessages(read_csv(master_key_path, col_types = cols()))

required_mk <- c(
  "GP2ID", "study", "wgs", "wgs_prune_reason", "wgs_label",
  "study_type", "baseline_GP2_phenotype", "biological_sex_for_qc",
  "age_at_sample_collection", "age_of_onset", "age_at_diagnosis",
  "age_at_death", "age_at_last_follow_up"
)
missing_mk <- setdiff(required_mk, names(mk))
if (length(missing_mk) > 0) {
  stop("Missing required columns in r10 master key: ", paste(missing_mk, collapse = ", "))
}

ext_clin <- suppressMessages(read_csv(ext_clin_path, col_types = cols()))

required_ext <- c(
  "GP2ID",
  "age_at_baseline",
  "age_at_onset",
  "age_at_first_motor_symptom",
  "age_at_diagnosis",
  "age_at_death",
  "age_at_living_status_censored",
  "age_at_dementia_status_censored",
  "age_brain_surgery"
)
missing_ext <- setdiff(required_ext, names(ext_clin))
if (length(missing_ext) > 0) {
  stop("Missing required columns in extended clinical data: ", paste(missing_ext, collapse = ", "))
}

# ------------------------------------------------------------------
# 3) Prepare join tables (with deterministic dedup for extended clinical)
# ------------------------------------------------------------------

mk_sub <- mk %>%
  transmute(
    GP2ID = as.character(GP2ID),
    study,
    wgs,
    wgs_prune_reason,
    wgs_label,
    study_type,
    baseline_GP2_phenotype,
    biological_sex_for_qc,
    age_at_sample_collection,
    age_of_onset,
    age_at_diagnosis,
    age_at_death,
    age_at_last_follow_up
  )

dup_mk <- mk_sub$GP2ID[duplicated(mk_sub$GP2ID) & !is.na(mk_sub$GP2ID)]
if (length(dup_mk) > 0) {
  stop(
    "Duplicate GP2ID values detected in master key (join would be ambiguous). ",
    "Examples: ", paste(unique(dup_mk)[seq_len(min(10, length(unique(dup_mk))))], collapse = ", ")
  )
}

# Extended clinical subset; rename overlaps so master key names remain unchanged
ext_sub_raw <- ext_clin %>%
  transmute(
    GP2ID = na_if(trimws(as.character(GP2ID)), ""),
    age_at_baseline,
    age_at_onset,
    age_at_first_motor_symptom,
    ext_age_at_diagnosis = age_at_diagnosis,
    ext_age_at_death     = age_at_death,
    age_at_living_status_censored,
    age_at_dementia_status_censored,
    age_brain_surgery
  )

cols_to_check <- setdiff(names(ext_sub_raw), "GP2ID")

# QC: duplicate GP2IDs (counts)
ext_dup_counts <- ext_sub_raw %>%
  filter(!is.na(GP2ID)) %>%
  count(GP2ID, name = "n_rows") %>%
  filter(n_rows > 1) %>%
  arrange(desc(n_rows), GP2ID)

if (nrow(ext_dup_counts) > 0) {
  cat("NOTE: Extended clinical data has duplicated GP2IDs.\n")
  cat("  Duplicated GP2IDs:", nrow(ext_dup_counts), "\n")
  cat("  Writing duplicate counts to:", qc_ext_dups, "\n\n")
  write_tsv(ext_dup_counts, qc_ext_dups)
}

# QC: conflicts (any field has >1 distinct non-NA value within a GP2ID)
# Use explicit column list to avoid dplyr grouped-across selection edge cases.
ext_conflicts <- ext_sub_raw %>%
  filter(!is.na(GP2ID)) %>%
  group_by(GP2ID) %>%
  summarise(
    across(all_of(cols_to_check), ~ n_distinct(na.omit(.x))),
    .groups = "drop"
  ) %>%
  # Any column count > 1 indicates conflicting values
  filter(if_any(all_of(cols_to_check), ~ .x > 1))

if (nrow(ext_conflicts) > 0) {
  cat("NOTE: Extended clinical data has GP2IDs with conflicting values across rows.\n")
  cat("  Conflicted GP2IDs:", nrow(ext_conflicts), "\n")
  cat("  Writing conflict summary to:", qc_ext_conf, "\n\n")
  write_tsv(ext_conflicts, qc_ext_conf)
}

# Deduplicate extended clinical to 1 row per GP2ID:
# Choose the row with the most non-missing values across selected fields.
ext_sub <- ext_sub_raw %>%
  filter(!is.na(GP2ID)) %>%
  mutate(.non_na = rowSums(!is.na(select(., all_of(cols_to_check))))) %>%
  group_by(GP2ID) %>%
  arrange(desc(.non_na)) %>%
  slice(1) %>%
  ungroup() %>%
  select(-.non_na)

# ------------------------------------------------------------------
# 4) Join core ped with master key + extended clinical
# ------------------------------------------------------------------

ped_annot <- core_ped %>%
  left_join(mk_sub,  by = c("id" = "GP2ID")) %>%
  left_join(ext_sub, by = c("id" = "GP2ID"))

# ------------------------------------------------------------------
# 5) Add convenience flags and coverage checks
# ------------------------------------------------------------------

ped_annot <- ped_annot %>%
  mutate(
    wgs_prune_reason_clean = trimws(coalesce(as.character(wgs_prune_reason), "")),
    has_wgs = case_when(
      is.na(wgs) ~ FALSE,
      as.character(wgs) == "1" & wgs_prune_reason_clean == "" ~ TRUE,
      TRUE ~ FALSE
    ),
    ancestry_label = wgs_label,
    pheno_from_master = case_when(
      baseline_GP2_phenotype %in% c("PD", "Affected_PD") ~ "PD",
      baseline_GP2_phenotype %in% c("Control", "Unaffected") ~ "Control",
      TRUE ~ "Other_or_missing"
    )
  ) %>%
  select(-wgs_prune_reason_clean)

n_no_master <- sum(is.na(ped_annot$study))
cat("Core ped individuals not found in r10 master key:", n_no_master, "\n")

n_with_wgs <- sum(ped_annot$has_wgs, na.rm = TRUE)
cat("Core ped individuals with WGS usable (wgs == 1 & no prune reason):", n_with_wgs, "\n")

n_with_ext_dx <- sum(!is.na(ped_annot$ext_age_at_diagnosis))
cat("Core ped individuals with extended age_at_diagnosis present:", n_with_ext_dx, "\n\n")

# ------------------------------------------------------------------
# 6) Write annotated core ped
# ------------------------------------------------------------------

write_tsv(ped_annot, out_ped_annot)

cat("Wrote annotated core r10 ped with master key and extended ages to:\n  ",
    out_ped_annot, "\n", sep = "")

