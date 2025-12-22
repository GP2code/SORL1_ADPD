#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ------------------------------------------------------------------
# Purpose:
#   Mark informative families and informative individuals in the merged
#   master pedigree (core + kinship).
#
# Informative family definition (as implemented here):
#   - >= 1 WGS PD individual in the family
#   - >= 2 WGS individuals total in the family
#
# Informative individual definition (as implemented here):
#   - In an informative family AND (has_wgs OR is_parent_of_WGS)
#
# Output location:
#   By default, writes outputs next to this script (script-local).
#   Override with OUT_DIR if needed.
#
# Env overrides:
#   OUT_DIR
#   MASTER_PED_PATH
#   MASTER_KEY_PATH
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

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------

master_ped_path <- Sys.getenv(
  "MASTER_PED_PATH",
  unset = file.path(out_dir, "r10_master_ped_all_families.tsv")
)

master_key_path <- Sys.getenv(
  "MASTER_KEY_PATH",
  unset = "/home/jupyter/workspace/gp2_tier2_eu_release10/clinical_data/master_key_release10_final_vwb.csv"
)

if (!file.exists(master_ped_path)) {
  stop(
    "Master PED not found: ", master_ped_path, "\n",
    "Run the merge script first in the same folder, or set MASTER_PED_PATH."
  )
}
if (!file.exists(master_key_path)) stop("Master key not found: ", master_key_path)

out_ped_inform <- file.path(out_dir, "r10_master_ped_informative_families.tsv")
out_fam_summary <- file.path(out_dir, "qc_informative_family_summary.tsv")

cat("Output directory (script-local by default):\n  ", out_dir, "\n\n", sep = "")

# ------------------------------------------------------------------
# 1) Read master ped (core + kinship) and r10 master key
# ------------------------------------------------------------------

ped <- read_tsv(
  master_ped_path,
  col_types = cols(
    family_id   = col_character(),
    id          = col_character(),
    parental_id = col_character(),
    maternal_id = col_character(),
    sex         = col_double(),
    phenotype   = col_double(),
    ancestry    = col_character(),
    source      = col_character()
  )
)

required_ped <- c("family_id", "id", "parental_id", "maternal_id", "sex", "phenotype")
missing_ped <- setdiff(required_ped, names(ped))
if (length(missing_ped) > 0) {
  stop("Master PED missing required columns: ", paste(missing_ped, collapse = ", "))
}

if (any(is.na(ped$family_id) | ped$family_id == "")) stop("Master PED has missing/blank family_id.")
if (any(is.na(ped$id) | ped$id == "")) stop("Master PED has missing/blank id.")

cat("r10 master ped (all families):\n")
cat("  Individuals:", nrow(ped), "\n")
cat("  Families   :", n_distinct(ped$family_id), "\n\n")

mk <- read_csv(master_key_path, col_types = cols())

required_mk <- c("GP2ID", "wgs", "wgs_prune_reason", "wgs_label", "baseline_GP2_phenotype")
missing_mk <- setdiff(required_mk, names(mk))
if (length(missing_mk) > 0) {
  stop("Missing required columns in r10 master key: ", paste(missing_mk, collapse = ", "))
}

mk_sub <- mk %>%
  transmute(
    GP2ID = as.character(GP2ID),
    wgs = as.character(wgs),
    wgs_prune_reason = as.character(wgs_prune_reason),
    wgs_label,
    baseline_GP2_phenotype
  )

dup_mk <- mk_sub$GP2ID[duplicated(mk_sub$GP2ID) & !is.na(mk_sub$GP2ID)]
if (length(dup_mk) > 0) {
  stop(
    "Duplicate GP2ID values detected in master key (join would be ambiguous). ",
    "Examples: ", paste(unique(dup_mk)[seq_len(min(10, length(unique(dup_mk))))], collapse = ", ")
  )
}

# ------------------------------------------------------------------
# 2) Join WGS info from master key onto ped + compute has_wgs (robust)
# ------------------------------------------------------------------

ped_annot <- ped %>%
  left_join(mk_sub, by = c("id" = "GP2ID")) %>%
  mutate(
    prune_clean = trimws(coalesce(wgs_prune_reason, "")),
    has_wgs = case_when(
      is.na(wgs) ~ FALSE,
      wgs == "1" & prune_clean == "" ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  select(-prune_clean)

cat("Individuals with usable WGS in master ped:", sum(ped_annot$has_wgs, na.rm = TRUE), "\n\n")

# ------------------------------------------------------------------
# 3) Define PD/control status from PED phenotype codes
#    phenotype: 2 = PD, 1 = control, 0/other = other/unknown
# ------------------------------------------------------------------

ped_annot <- ped_annot %>%
  mutate(
    phenotype_i = suppressWarnings(as.integer(phenotype)),
    is_PD       = phenotype_i == 2L,
    is_ctrl     = phenotype_i == 1L,
    is_WGS_PD   = has_wgs & is_PD,
    is_WGS_ctrl = has_wgs & is_ctrl
  ) %>%
  select(-phenotype_i)

# ------------------------------------------------------------------
# 4) Summarize per family and mark informative families
#    Informative family: >=1 WGS PD and >=2 WGS total
# ------------------------------------------------------------------

fam_summary <- ped_annot %>%
  group_by(family_id) %>%
  summarise(
    n_members   = n(),
    n_WGS_PD    = sum(is_WGS_PD, na.rm = TRUE),
    n_WGS_ctrl  = sum(is_WGS_ctrl, na.rm = TRUE),
    n_WGS_total = sum(has_wgs, na.rm = TRUE),
    informative_family = (n_WGS_PD >= 1 & n_WGS_total >= 2),
    .groups = "drop"
  ) %>%
  arrange(desc(informative_family), desc(n_WGS_PD), desc(n_WGS_total), family_id)

cat("Total families:", nrow(fam_summary), "\n")
cat("Informative families (>=1 WGS PD & >=2 WGS total):",
    sum(fam_summary$informative_family), "\n\n")

write_tsv(fam_summary, out_fam_summary)
cat("Wrote family-level QC summary to:\n  ", out_fam_summary, "\n\n", sep = "")

# ------------------------------------------------------------------
# 5) Mark informative individuals
#    - informative_family = TRUE at family level
#    - is_informative = TRUE if (has_wgs OR is_parent_of_WGS) in an informative family
# ------------------------------------------------------------------

ped2 <- ped_annot %>%
  left_join(fam_summary %>% select(family_id, informative_family), by = "family_id") %>%
  mutate(informative_family = if_else(is.na(informative_family), FALSE, informative_family))

# Identify parents (present in PED) of WGS individuals within informative families
wgs_inds <- ped2 %>%
  filter(informative_family, has_wgs) %>%
  select(family_id, parental_id, maternal_id)

parents_to_flag <- ped2 %>%
  filter(informative_family) %>%
  semi_join(wgs_inds, by = c("family_id", "id" = "parental_id")) %>%
  bind_rows(
    ped2 %>%
      filter(informative_family) %>%
      semi_join(wgs_inds, by = c("family_id", "id" = "maternal_id"))
  ) %>%
  distinct(family_id, id) %>%
  mutate(is_parent_of_WGS = TRUE)

ped3 <- ped2 %>%
  left_join(parents_to_flag, by = c("family_id", "id")) %>%
  mutate(
    is_parent_of_WGS = if_else(is.na(is_parent_of_WGS), FALSE, is_parent_of_WGS),
    is_informative   = informative_family & (has_wgs | is_parent_of_WGS)
  ) %>%
  arrange(family_id, id)

# ------------------------------------------------------------------
# 6) Counts + write output
# ------------------------------------------------------------------

cat("Total informative individuals (is_informative == TRUE):",
    sum(ped3$is_informative, na.rm = TRUE), "\n")
cat("  of which WGS:", sum(ped3$is_informative & ped3$has_wgs, na.rm = TRUE), "\n")
cat("  of which non-WGS parents:",
    sum(ped3$is_informative & !ped3$has_wgs, na.rm = TRUE), "\n\n")

write_tsv(ped3, out_ped_inform)
cat("Wrote r10 master ped with informative flags to:\n  ", out_ped_inform, "\n", sep = "")

