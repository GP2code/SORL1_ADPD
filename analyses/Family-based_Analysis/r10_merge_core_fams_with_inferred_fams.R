#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ------------------------------------------------------------------
# Purpose:
#   Merge the annotated "core" pedigree with the kinship-inferred pedigree
#   into a unified master PED-like table (all families).
#
# Output location:
#   By default, writes outputs next to this script (script-local).
#   Override with OUT_DIR if needed.
#
# Env overrides:
#   OUT_DIR
#   CORE_PED_ANNOT_PATH
#   KINSHIP_PED_PATH
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

# Default inputs: script-local outputs from Scripts 1 and 2
core_ped_annot_path <- Sys.getenv(
  "CORE_PED_ANNOT_PATH",
  unset = file.path(out_dir, "r10_core_ped_with_master_ages.tsv")
)

kinship_ped_path <- Sys.getenv(
  "KINSHIP_PED_PATH",
  unset = file.path(out_dir, "r10_kinship_inferred_families_ped.tsv")
)

if (!file.exists(core_ped_annot_path)) {
  stop(
    "Core annotated PED not found: ", core_ped_annot_path, "\n",
    "Run Script 1 in the same folder first, or set CORE_PED_ANNOT_PATH."
  )
}
if (!file.exists(kinship_ped_path)) {
  stop(
    "Kinship-inferred PED not found: ", kinship_ped_path, "\n",
    "Run Script 2 in the same folder first, or set KINSHIP_PED_PATH."
  )
}

out_master_ped <- file.path(out_dir, "r10_master_ped_all_families.tsv")
out_family_sizes <- file.path(out_dir, "qc_master_ped_family_sizes.tsv")
out_src_by_anc   <- file.path(out_dir, "qc_master_ped_source_by_ancestry.tsv")

cat("Output directory (script-local by default):\n  ", out_dir, "\n\n", sep = "")

# ------------------------------------------------------------------
# 1) Read inputs
# ------------------------------------------------------------------

core <- read_tsv(core_ped_annot_path, col_types = cols())
kin  <- read_tsv(kinship_ped_path,    col_types = cols())

required_core_cols <- c("family_id", "id", "parental_id", "maternal_id", "sex", "phenotype")
missing_core <- setdiff(required_core_cols, names(core))
if (length(missing_core) > 0) {
  stop("Core annotated ped is missing required columns: ", paste(missing_core, collapse = ", "))
}

required_kin_cols <- c("family_id", "id", "parental_id", "maternal_id", "sex", "phenotype", "ancestry", "source")
missing_kin <- setdiff(required_kin_cols, names(kin))
if (length(missing_kin) > 0) {
  stop("Kinship ped is missing required columns: ", paste(missing_kin, collapse = ", "))
}

cat("Core ped annotated:\n")
cat("  Individuals:", nrow(core), "\n")
cat("  Families   :", n_distinct(core$family_id), "\n\n")

cat("Kinship-inferred ped:\n")
cat("  Individuals:", nrow(kin), "\n")
cat("  Families   :", n_distinct(kin$family_id), "\n\n")

# ------------------------------------------------------------------
# 2) Determine ancestry column from core annotated PED
# ------------------------------------------------------------------

core_anc_col <- NULL
if ("ancestry" %in% names(core)) {
  core_anc_col <- "ancestry"
} else if ("ancestry_label" %in% names(core)) {
  core_anc_col <- "ancestry_label"
} else if ("wgs_label" %in% names(core)) {
  core_anc_col <- "wgs_label"
}

if (is.null(core_anc_col)) {
  stop("Core annotated ped does not contain an ancestry column (expected one of: ancestry, ancestry_label, wgs_label).")
}

# ------------------------------------------------------------------
# 3) Validate uniqueness / overlap constraints
# ------------------------------------------------------------------

dup_core <- core %>%
  count(family_id, id, name = "n") %>%
  filter(n > 1)
if (nrow(dup_core) > 0) {
  stop("Duplicate (family_id, id) rows found in core ped annotated (unexpected).")
}

dup_kin <- kin %>%
  count(family_id, id, name = "n") %>%
  filter(n > 1)
if (nrow(dup_kin) > 0) {
  stop("Duplicate (family_id, id) rows found in kinship ped (unexpected).")
}

# IDs should not overlap (Script 2 intentionally excluded core IDs)
overlap_ids <- intersect(unique(core$id), unique(kin$id))
if (length(overlap_ids) > 0) {
  stop(
    "Overlap detected between core IDs and kinship IDs (unexpected). ",
    "Examples: ", paste(overlap_ids[seq_len(min(10, length(overlap_ids)))], collapse = ", ")
  )
}

# ------------------------------------------------------------------
# 4) Prepare unified schema
# ------------------------------------------------------------------

fix_parent <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "0"
  x
}

fix_int_code <- function(x, default = 0L) {
  y <- suppressWarnings(as.integer(x))
  y[is.na(y)] <- default
  y
}

core_out <- core %>%
  transmute(
    family_id   = as.character(family_id),
    id          = as.character(id),
    parental_id = fix_parent(parental_id),
    maternal_id = fix_parent(maternal_id),
    sex         = fix_int_code(sex, default = 0L),
    phenotype   = fix_int_code(phenotype, default = 0L),
    ancestry    = as.character(.data[[core_anc_col]]),
    source      = "core"
  )

kin_out <- kin %>%
  transmute(
    family_id   = as.character(family_id),
    id          = as.character(id),
    parental_id = fix_parent(parental_id),
    maternal_id = fix_parent(maternal_id),
    sex         = fix_int_code(sex, default = 0L),
    phenotype   = fix_int_code(phenotype, default = 0L),
    ancestry    = as.character(ancestry),
    source      = as.character(source)
  )

# Sanity: ensure no missing family_id/id
if (any(is.na(core_out$family_id) | core_out$family_id == "") ||
    any(is.na(core_out$id) | core_out$id == "")) {
  stop("Core output contains missing family_id or id after standardization.")
}
if (any(is.na(kin_out$family_id) | kin_out$family_id == "") ||
    any(is.na(kin_out$id) | kin_out$id == "")) {
  stop("Kinship output contains missing family_id or id after standardization.")
}

# ------------------------------------------------------------------
# 5) Combine into master ped (stable row order)
# ------------------------------------------------------------------

master_ped <- bind_rows(core_out, kin_out) %>%
  arrange(family_id, id)

# Final check: each id must belong to exactly 1 family
multi_fam_ids <- master_ped %>%
  distinct(id, family_id) %>%
  count(id, name = "n_fams") %>%
  filter(n_fams > 1)

if (nrow(multi_fam_ids) > 0) {
  stop("Some IDs appear in multiple families in the combined master ped (unexpected).")
}

cat("Combined master ped:\n")
cat("  Individuals:", nrow(master_ped), "\n")
cat("  Families   :", n_distinct(master_ped$family_id), "\n\n")

cat("  Core individuals   :", sum(master_ped$source == "core"), "\n")
cat("  Kinship individuals:", sum(master_ped$source != "core"), "\n\n")

# ------------------------------------------------------------------
# 6) Write master ped and QC summaries
# ------------------------------------------------------------------

write_tsv(master_ped, out_master_ped)
cat("Wrote r10 master ped (all families) to:\n  ", out_master_ped, "\n", sep = "")

family_sizes <- master_ped %>%
  count(family_id, name = "n_members") %>%
  count(n_members, name = "n_families") %>%
  arrange(n_members)

src_by_anc <- master_ped %>%
  count(ancestry, source) %>%
  arrange(ancestry, source)

write_tsv(family_sizes, out_family_sizes)
write_tsv(src_by_anc,   out_src_by_anc)

cat("\nFamily size distribution written to:\n  ", out_family_sizes, "\n", sep = "")
cat("Source breakdown by ancestry written to:\n  ", out_src_by_anc, "\n\n", sep = "")

cat("Family size distribution (n_members per family):\n")
print(family_sizes)

cat("\nSource breakdown by ancestry:\n")
print(src_by_anc)

