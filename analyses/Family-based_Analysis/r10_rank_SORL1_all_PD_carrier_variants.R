#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ------------------------------------------------------------------
# Purpose:
#   Rank SORL1 variants by the number of families in which:
#     - at least one PD is genotyped (n_PD > 0)
#     - all genotyped PD are carriers (n_noncarrier_PD == 0)
#     - and at least one PD carrier exists (n_carrier_PD > 0)
#
# NEW stricter criterion requested:
#   Exclude variants that also "appear in other families" where not all PD
#   are carriers. Operationally:
#     - define fam_PD_violation = (n_PD>0 & variant_present_in_family & n_noncarrier_PD>0)
#     - keep only variants with n_families_PD_violation == 0
#
# Output location:
#   By default, reads/writes next to this script (script-local).
#
# Env overrides:
#   OUT_DIR
#   BY_FAMILY_PATH
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

by_family_path <- Sys.getenv(
  "BY_FAMILY_PATH",
  unset = file.path(out_dir, "r10_SORL1_segregation_by_family.tsv")
)

out_ranked_strict <- file.path(out_dir, "r10_SORL1_variants_all_PD_carriers_ranked.tsv")
out_fam_hits      <- file.path(out_dir, "qc_SORL1_families_all_PD_carriers.tsv")
out_fam_viol      <- file.path(out_dir, "qc_SORL1_families_PD_violation.tsv")

if (!file.exists(by_family_path)) {
  stop(
    "Per-family segregation table not found: ", by_family_path, "\n",
    "Run r10_summarize_SORL1_segregation_by_ancestry.R first (in the same folder), or set BY_FAMILY_PATH."
  )
}

cat("OUT_DIR        :", out_dir, "\n")
cat("BY_FAMILY_PATH :", by_family_path, "\n\n")

# ------------------------------------------------------------------
# 1) Read per-family segregation file
# ------------------------------------------------------------------

seg_fam <- read_tsv(by_family_path, col_types = cols())

required_cols <- c(
  "ancestry", "variant_id", "family_id",
  "n_PD", "n_ctrl", "n_other",
  "n_carrier_PD", "n_carrier_ctrl", "n_carrier_other",
  "n_noncarrier_PD", "n_noncarrier_ctrl", "n_noncarrier_other"
)

missing_cols <- setdiff(required_cols, names(seg_fam))
if (length(missing_cols) > 0) {
  stop("Missing required columns in r10_SORL1_segregation_by_family.tsv: ",
       paste(missing_cols, collapse = ", "))
}

cat("Rows in per-family segregation table:", nrow(seg_fam), "\n\n")

# Ensure numeric columns are numeric
num_cols <- setdiff(required_cols, c("ancestry", "variant_id", "family_id"))
seg_fam <- seg_fam %>%
  mutate(across(all_of(num_cols), ~ suppressWarnings(as.integer(.x))))

# ------------------------------------------------------------------
# 2) Family-level flags
# ------------------------------------------------------------------

seg_flags <- seg_fam %>%
  mutate(
    # Any carrier in the family for this variant (PD/ctrl/other)
    fam_has_any_carrier = (n_carrier_PD + n_carrier_ctrl + n_carrier_other) > 0,

    # All genotyped PD are carriers (and at least one PD carrier)
    fam_all_PD_carriers = (n_PD > 0 & n_carrier_PD > 0 & n_noncarrier_PD == 0),

    # PD-violation family: variant present in family, but not all genotyped PD are carriers
    fam_PD_violation = (n_PD > 0 & fam_has_any_carrier & n_noncarrier_PD > 0),

    fam_any_ctrl_carrier = (n_carrier_ctrl > 0)
  )

# QC: families meeting the "all PD carriers" criterion
fam_hits <- seg_flags %>%
  filter(fam_all_PD_carriers) %>%
  arrange(ancestry, variant_id, family_id)

write_tsv(fam_hits, out_fam_hits)

cat("Families meeting all-PD-carriers criterion (rows):", nrow(fam_hits), "\n")
cat("Wrote family-level hit table to:\n  ", out_fam_hits, "\n\n", sep = "")

# QC: families that violate strict segregation (variant present, PD noncarriers exist)
fam_viol <- seg_flags %>%
  filter(fam_PD_violation) %>%
  arrange(ancestry, variant_id, family_id)

write_tsv(fam_viol, out_fam_viol)

cat("Families with PD-violation (variant present but some PD noncarriers) rows:", nrow(fam_viol), "\n")
cat("Wrote PD-violation family table to:\n  ", out_fam_viol, "\n\n", sep = "")

# ------------------------------------------------------------------
# 3) Variant-level ranking summary + strict exclusion of PD-violation variants
# ------------------------------------------------------------------

var_rank <- seg_flags %>%
  group_by(ancestry, variant_id) %>%
  summarise(
    n_families_total           = n_distinct(family_id),

    n_families_with_any_PD     = sum(n_PD > 0, na.rm = TRUE),
    n_families_with_PD_carriers= sum(n_carrier_PD > 0, na.rm = TRUE),

    n_families_all_PD_carriers = sum(fam_all_PD_carriers, na.rm = TRUE),
    n_families_PD_violation    = sum(fam_PD_violation, na.rm = TRUE),

    n_families_all_PD_carriers_no_ctrl_carriers =
      sum(fam_all_PD_carriers & !fam_any_ctrl_carrier, na.rm = TRUE),

    # Totals restricted to "all-PD-carrier families"
    total_PD_carriers_all_PDfam    = sum(n_carrier_PD[fam_all_PD_carriers], na.rm = TRUE),
    total_ctrl_carriers_all_PDfam  = sum(n_carrier_ctrl[fam_all_PD_carriers], na.rm = TRUE),
    total_other_carriers_all_PDfam = sum(n_carrier_other[fam_all_PD_carriers], na.rm = TRUE),
    total_PD_nonmissing_all_PDfam  = sum(n_PD[fam_all_PD_carriers], na.rm = TRUE),
    total_ctrl_nonmissing_all_PDfam= sum(n_ctrl[fam_all_PD_carriers], na.rm = TRUE),
    .groups = "drop"
  )

# Apply strict filter:
#   - must have >=1 all-PD-carrier family
#   - must have 0 PD-violation families
var_strict <- var_rank %>%
  filter(n_families_all_PD_carriers >= 1, n_families_PD_violation == 0) %>%
  arrange(
    ancestry,
    desc(n_families_all_PD_carriers),
    desc(total_PD_carriers_all_PDfam),
    total_ctrl_carriers_all_PDfam,   # ascending default
    variant_id
  )

write_tsv(var_strict, out_ranked_strict)

cat("Wrote STRICT ranked variants to:\n  ", out_ranked_strict, "\n", sep = "")
cat("Number of variants passing strict criteria:", nrow(var_strict), "\n")

n_candidates_relaxed <- sum(var_rank$n_families_all_PD_carriers >= 1, na.rm = TRUE)
n_removed_by_violation <- sum((var_rank$n_families_all_PD_carriers >= 1) & (var_rank$n_families_PD_violation > 0), na.rm = TRUE)

cat("Variants with >=1 all-PD-carrier family (relaxed):", n_candidates_relaxed, "\n")
cat("Of those, removed due to PD-violation families:", n_removed_by_violation, "\n")

