#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(igraph)
})

# ------------------------------------------------------------------
# Purpose:
#   Build kinship-inferred families for monogenic WGS samples not in the
#   "core" pedigree, using ancestry-specific PLINK2 .related files.
#
# Output location:
#   By default, writes outputs next to this script (script-local).
#   Override with OUT_DIR if needed.
#
# Env overrides:
#   OUT_DIR
#   CORE_PED_ANNOT_PATH
#   MASTER_KEY_PATH
#   RELATED_BASE
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

core_ped_annot_path <- Sys.getenv(
  "CORE_PED_ANNOT_PATH",
  unset = file.path(out_dir, "r10_core_ped_with_master_ages.tsv")
)

master_key_path <- Sys.getenv(
  "MASTER_KEY_PATH",
  unset = "/home/jupyter/workspace/gp2_tier2_eu_release10/clinical_data/master_key_release10_final_vwb.csv"
)

related_base <- Sys.getenv(
  "RELATED_BASE",
  unset = "/home/jupyter/workspace/gp2_tier2_eu_release10/wgs/deepvariant_joint_calling/related_samples"
)

if (!file.exists(core_ped_annot_path)) {
  stop(
    "Annotated core PED not found: ", core_ped_annot_path, "\n",
    "Run Script 1 in the same folder first, or set CORE_PED_ANNOT_PATH."
  )
}
if (!file.exists(master_key_path)) stop("Master key not found: ", master_key_path)
if (!dir.exists(related_base)) stop("RELATED_BASE directory not found: ", related_base)

out_kinship_ped <- file.path(out_dir, "r10_kinship_inferred_families_ped.tsv")
out_singletons  <- file.path(out_dir, "qc_monogenic_noncore_singletons.tsv")

cat("Output directory (script-local by default):\n  ", out_dir, "\n\n", sep = "")

# ------------------------------------------------------------------
# 1) Read core ped (annotated) and master key
# ------------------------------------------------------------------

core_ped <- read_tsv(core_ped_annot_path, col_types = cols())

if (!all(c("family_id", "id") %in% names(core_ped))) {
  stop("core_ped missing required columns: family_id, id")
}

cat("Core r10 ped annotated:\n")
cat("  Individuals:", nrow(core_ped), "\n")
cat("  Families   :", n_distinct(core_ped$family_id), "\n\n")

core_ids <- unique(core_ped$id)

mk <- read_csv(master_key_path, col_types = cols())

required_mk <- c(
  "GP2ID", "study", "wgs", "wgs_prune_reason",
  "wgs_label", "study_type",
  "baseline_GP2_phenotype", "biological_sex_for_qc",
  "age_at_sample_collection"
)
missing_mk <- setdiff(required_mk, names(mk))
if (length(missing_mk) > 0) {
  stop("Missing required columns in r10 master key: ", paste(missing_mk, collapse = ", "))
}

# ------------------------------------------------------------------
# 2) Define monogenic WGS individuals (R10) excluding core
# ------------------------------------------------------------------

mk2 <- mk %>%
  mutate(
    GP2ID = as.character(GP2ID),
    wgs_chr = as.character(wgs),
    prune_clean = trimws(coalesce(as.character(wgs_prune_reason), "")),
    study_type_clean = trimws(coalesce(as.character(study_type), ""))
  )

mono_wgs <- mk2 %>%
  filter(
    wgs_chr == "1",
    prune_clean == "",
    study_type_clean == "Monogenic"
  ) %>%
  mutate(is_core = GP2ID %in% core_ids)

cat("Total Monogenic WGS samples (released, QC-passed):", nrow(mono_wgs), "\n")
cat("  Of these, in core ped:", sum(mono_wgs$is_core), "\n")
cat("  To be used for kinship-inferred families:", sum(!mono_wgs$is_core), "\n\n")

mono_wgs_new <- mono_wgs %>%
  filter(!is_core)

# ------------------------------------------------------------------
# 3) Helper: phenotype and sex coding for PED
# ------------------------------------------------------------------

ped_pheno_code <- function(pheno) {
  dplyr::case_when(
    pheno %in% c("PD", "Affected_PD") ~ 2L,
    pheno %in% c("Control", "Unaffected") ~ 1L,
    TRUE ~ 0L
  )
}

ped_sex_code <- function(sex_str) {
  dplyr::case_when(
    sex_str == "Male"   ~ 1L,
    sex_str == "Female" ~ 2L,
    TRUE                ~ 0L
  )
}

# ------------------------------------------------------------------
# 4) Build kinship clusters and convert to PED per ancestry (deterministic)
# ------------------------------------------------------------------

ancestries <- mono_wgs_new %>%
  filter(!is.na(wgs_label), trimws(wgs_label) != "") %>%
  pull(wgs_label) %>%
  unique() %>%
  sort()

cat("Monogenic WGS non-core ancestries:", paste(ancestries, collapse = ", "), "\n\n")

kinship_ped_list <- list()
singleton_list   <- list()

for (anc in ancestries) {
  cat("=== Processing ancestry:", anc, "===\n")
  rel_file <- file.path(related_base, anc, paste0(anc, "_release10.related"))

  mono_anc <- mono_wgs_new %>%
    filter(wgs_label == anc)

  cat("  Monogenic WGS non-core samples in", anc, ":", nrow(mono_anc), "\n")

  if (nrow(mono_anc) == 0) {
    cat("  No non-core monogenic samples for", anc, " – skipping.\n\n")
    next
  }

  # IDs in this ancestry set
  anc_ids <- sort(unique(mono_anc$GP2ID))

  if (!file.exists(rel_file)) {
    cat("  No related file found at:", rel_file, "– skipping kinship inference for this ancestry.\n")
    # Still record them as "no related file"
    singleton_list[[length(singleton_list) + 1L]] <- mono_anc %>%
      transmute(
        ancestry = anc,
        GP2ID,
        baseline_GP2_phenotype,
        reason = "missing_related_file"
      )
    cat("\n")
    next
  }

  rel <- read_delim(
    rel_file,
    delim = ",",
    comment = "#",
    col_names = c("FID1","IID1","FID2","IID2","NSNP","HETHET","IBS0","KINSHIP","REL"),
    col_types = cols(
      FID1 = col_character(),
      IID1 = col_character(),
      FID2 = col_character(),
      IID2 = col_character(),
      NSNP = col_double(),
      HETHET = col_double(),
      IBS0 = col_double(),
      KINSHIP = col_double(),
      REL = col_character()
    )
  )

  rel_sub <- rel %>%
    filter(IID1 %in% anc_ids & IID2 %in% anc_ids) %>%
    filter(IID1 != IID2) %>%
    select(IID1, IID2, REL, KINSHIP) %>%
    distinct()

  cat("  Related pairs within monogenic non-core", anc, ":", nrow(rel_sub), "\n")

  # Build graph INCLUDING all anc_ids as vertices (so singletons are tracked)
  g <- graph_from_data_frame(
    rel_sub[, c("IID1", "IID2")],
    directed = FALSE,
    vertices = data.frame(name = anc_ids)
  )
  comps <- components(g)
  membership <- comps$membership
  vertex_ids <- names(membership)  # should include all anc_ids

  # Split into clusters by component
  clusters <- split(vertex_ids, membership)

  # Identify singletons (clusters of size 1)
  singleton_ids <- sort(unlist(clusters[vapply(clusters, length, integer(1)) == 1], use.names = FALSE))
  if (length(singleton_ids) > 0) {
    singleton_list[[length(singleton_list) + 1L]] <- mono_anc %>%
      filter(GP2ID %in% singleton_ids) %>%
      transmute(
        ancestry = anc,
        GP2ID,
        baseline_GP2_phenotype,
        reason = "no_kinship_edges"
      )
  }

  # Keep only clusters of size >=2 for new kinship families
  family_clusters <- clusters[vapply(clusters, length, integer(1)) >= 2]

  if (length(family_clusters) == 0) {
    cat("  No clusters (size>=2) for", anc, " – no new families inferred.\n\n")
    next
  }

  # Deterministic ordering:
  # sort members within each cluster; then sort clusters by their smallest member ID
  family_clusters_sorted <- lapply(family_clusters, function(v) sort(v))
  cluster_keys <- vapply(family_clusters_sorted, function(v) v[1], character(1))
  ord <- order(cluster_keys)
  family_clusters_sorted <- family_clusters_sorted[ord]

  mono_anc_df <- mono_anc %>%
    select(
      GP2ID, study, wgs_label, baseline_GP2_phenotype,
      biological_sex_for_qc, age_at_sample_collection
    )

  fam_counter <- 0L
  for (cluster_vertices in family_clusters_sorted) {
    fam_counter <- fam_counter + 1L
    fam_id <- sprintf("%s_REL_FAM%04d", anc, fam_counter)

    cluster_df <- mono_anc_df %>%
      filter(GP2ID %in% cluster_vertices)

    fam_ped <- cluster_df %>%
      transmute(
        family_id   = fam_id,
        id          = GP2ID,
        parental_id = "0",
        maternal_id = "0",
        sex         = ped_sex_code(biological_sex_for_qc),
        phenotype   = ped_pheno_code(baseline_GP2_phenotype),
        ancestry    = wgs_label,
        source      = "related_kinship_r10"
      )

    kinship_ped_list[[length(kinship_ped_list) + 1L]] <- fam_ped
  }

  cat("  New kinship-inferred families created for", anc, ":", fam_counter, "\n\n")
}

# ------------------------------------------------------------------
# 5) Write outputs
# ------------------------------------------------------------------

# Always write singleton QC file (even if empty) for transparency
singletons_df <- if (length(singleton_list) == 0) {
  tibble(
    ancestry = character(),
    GP2ID = character(),
    baseline_GP2_phenotype = character(),
    reason = character()
  )
} else {
  bind_rows(singleton_list)
}
write_tsv(singletons_df, out_singletons)
cat("Wrote singleton/QC list to:\n  ", out_singletons, "\n\n", sep = "")

# Always write kinship PED (even if empty) so downstream scripts can rely on it
kinship_ped <- if (length(kinship_ped_list) == 0) {
  tibble(
    family_id = character(),
    id = character(),
    parental_id = character(),
    maternal_id = character(),
    sex = integer(),
    phenotype = integer(),
    ancestry = character(),
    source = character()
  )
} else {
  bind_rows(kinship_ped_list)
}

cat("Total kinship-inferred families (all ancestries):", n_distinct(kinship_ped$family_id), "\n")
cat("Total individuals in kinship-inferred families:", nrow(kinship_ped), "\n\n")

write_tsv(kinship_ped, out_kinship_ped)
cat("Wrote kinship-inferred r10 families PED to:\n  ", out_kinship_ped, "\n", sep = "")

