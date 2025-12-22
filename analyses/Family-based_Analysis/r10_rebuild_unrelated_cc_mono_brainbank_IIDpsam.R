#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(igraph)
})

Sys.setenv(TZ = "UTC")
options(readr.show_col_types = FALSE)

# ------------------------------------------------------------------
# Script-local output by default (can override with OUT_DIR)
# ------------------------------------------------------------------
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 0L) return(getwd())
  script_path <- sub("^--file=", "", file_arg[1])
  dirname(normalizePath(script_path))
}

script_dir <- get_script_dir()

# ------------------------------------------------------------------
# Inputs (env overrides supported)
# ------------------------------------------------------------------
master_key_csv <- Sys.getenv(
  "MASTER_KEY_CSV",
  unset = "/home/jupyter/workspace/gp2_tier2_eu_release10/clinical_data/master_key_release10_final_vwb.csv"
)
related_base <- Sys.getenv(
  "RELATED_BASE",
  unset = "/home/jupyter/workspace/gp2_tier2_eu_release10/wgs/deepvariant_joint_calling/related_samples"
)
plink_base <- Sys.getenv(
  "PLINK_BASE",
  unset = "/home/jupyter/workspace/gp2_tier2_eu_release10/wgs/deepvariant_joint_calling/plink"
)

# Output base directory (script-local by default)
out_dir <- Sys.getenv("OUT_DIR", unset = file.path(script_dir, "cc_data"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------
# Outputs
# ------------------------------------------------------------------
out_pheno   <- file.path(out_dir, "r10_wgs_cc_mono_plus_brainbank_UNRELATED_phenotype.tsv")
out_keep    <- file.path(out_dir, "r10_wgs_cc_mono_plus_brainbank_UNRELATED_plink_keep.tsv")
out_iids    <- file.path(out_dir, "r10_wgs_cc_mono_plus_brainbank_UNRELATED_IID_only.txt")
out_counts  <- file.path(out_dir, "r10_wgs_cc_mono_plus_brainbank_UNRELATED_casecontrol_counts_by_ancestry.tsv")

qc_prune_summary <- file.path(out_dir, "qc_cc_unrelated_prune_summary_by_ancestry.tsv")
qc_selection_tbl <- file.path(out_dir, "qc_cc_unrelated_selection_table_prepsam.tsv")
qc_psam_used     <- file.path(out_dir, "qc_cc_psam_files_used_by_ancestry.tsv")

cat("OUT_DIR           :", out_dir, "\n")
cat("MASTER_KEY_CSV    :", master_key_csv, "\n")
cat("RELATED_BASE      :", related_base, "\n")
cat("PLINK_BASE        :", plink_base, "\n\n")

if (!file.exists(master_key_csv)) stop("Master key CSV not found: ", master_key_csv)
if (!dir.exists(related_base)) cat("NOTE: RELATED_BASE dir not found (script will treat all as unrelated if ancestry files missing).\n")
if (!dir.exists(plink_base)) stop("PLINK_BASE dir not found: ", plink_base)

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
clean_chr <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, '"', "")
  x <- str_trim(x)
  x[x == ""] <- NA_character_
  x[tolower(x) %in% c("na", "nan")] <- NA_character_
  x
}

# Robust PSAM reader: pick any chr*_ANC_release10.psam deterministically
read_psam <- function(anc) {
  # Prefer chr11 then chr1, else first in sorted glob
  candidates <- c(
    file.path(plink_base, anc, paste0("chr11_", anc, "_release10.psam")),
    file.path(plink_base, anc, paste0("chr1_",  anc, "_release10.psam"))
  )

  psam_path <- NA_character_
  for (p in candidates) {
    if (file.exists(p)) { psam_path <- p; break }
  }

  if (is.na(psam_path)) {
    glob <- Sys.glob(file.path(plink_base, anc, paste0("chr*_", anc, "_release10.psam")))
    if (length(glob) == 0) return(NULL)
    psam_path <- sort(glob)[1]
  }

  ps <- read_tsv(psam_path, col_types = cols(), progress = FALSE)
  names(ps)[1] <- sub("^#", "", names(ps)[1])
  if (!all(c("FID", "IID") %in% names(ps))) {
    stop("PSAM missing FID/IID: ", psam_path)
  }

  list(path = psam_path, data = ps %>% select(FID, IID))
}

choose_representative <- function(cluster_df) {
  # cluster_df must contain: GP2ID, PHENO_PLINK, age_onset_or_recruitment
  cases <- cluster_df %>% filter(PHENO_PLINK == 2L)
  ctrls <- cluster_df %>% filter(PHENO_PLINK == 1L)

  if (nrow(cases) > 0) {
    cases_nonNA <- cases %>% filter(!is.na(age_onset_or_recruitment))
    if (nrow(cases_nonNA) > 0) {
      return(cases_nonNA %>% arrange(age_onset_or_recruitment) %>% slice(1) %>% pull(GP2ID))
    }
    return(sample(cases$GP2ID, 1))
  }

  if (nrow(ctrls) > 0) {
    ctrls_nonNA <- ctrls %>% filter(!is.na(age_onset_or_recruitment))
    if (nrow(ctrls_nonNA) > 0) {
      return(ctrls_nonNA %>% arrange(desc(age_onset_or_recruitment)) %>% slice(1) %>% pull(GP2ID))
    }
    return(sample(ctrls$GP2ID, 1))
  }

  sample(cluster_df$GP2ID, 1)
}

infer_unrelated_selection <- function(dat_anc, rel_file) {
  ids <- dat_anc$GP2ID

  base_tbl <- dat_anc %>%
    select(GP2ID, PHENO_PLINK, PHENO, age_onset_or_recruitment) %>%
    mutate(
      in_graph = FALSE,
      cluster_id = NA_integer_,
      cluster_size = 1L,
      representative = GP2ID,
      KEEP_UNRELATED = TRUE,
      note = NA_character_
    )

  if (!file.exists(rel_file)) {
    base_tbl$note <- "no_related_file"
    return(list(
      keep_ids = ids,
      n_edges = 0L,
      n_clusters_gt1 = 0L,
      selection_tbl = base_tbl
    ))
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
    ),
    progress = FALSE
  )

  edges <- rel %>%
    filter(IID1 %in% ids, IID2 %in% ids) %>%
    select(IID1, IID2) %>%
    distinct()

  if (nrow(edges) == 0) {
    base_tbl$note <- "no_pairs_in_selected_set"
    return(list(
      keep_ids = ids,
      n_edges = 0L,
      n_clusters_gt1 = 0L,
      selection_tbl = base_tbl
    ))
  }

  g <- graph_from_data_frame(edges, directed = FALSE)
  comps <- components(g)

  membership <- comps$membership           # named by IID
  csize <- comps$csize                     # indexed by component id
  n_clusters_gt1 <- sum(csize > 1)

  # Determine representative per component
  rep_by_cluster <- rep(NA_character_, length(csize))
  for (comp_id in seq_along(csize)) {
    if (csize[comp_id] <= 1) next
    members <- names(membership)[membership == comp_id]
    cl_df <- dat_anc %>% filter(GP2ID %in% members)
    rep_by_cluster[comp_id] <- choose_representative(cl_df)
  }

  # Build selection table
  cid <- unname(membership[base_tbl$GP2ID])
  base_tbl <- base_tbl %>%
    mutate(
      in_graph = !is.na(cid),
      cluster_id = cid,
      cluster_size = if_else(is.na(cluster_id), 1L, as.integer(csize[cluster_id])),
      representative = if_else(
        is.na(cluster_id),
        GP2ID,
        rep_by_cluster[cluster_id]
      ),
      KEEP_UNRELATED = if_else(
        is.na(cluster_id),
        TRUE,
        GP2ID == representative
      ),
      note = if_else(is.na(cluster_id), "singleton_or_not_in_related_edges", "in_related_component")
    )

  keep_ids <- base_tbl %>% filter(KEEP_UNRELATED) %>% pull(GP2ID)

  list(
    keep_ids = keep_ids,
    n_edges = nrow(edges),
    n_clusters_gt1 = n_clusters_gt1,
    selection_tbl = base_tbl
  )
}

# ------------------------------------------------------------------
# 1) Read master key
# ------------------------------------------------------------------
mk <- read_csv(master_key_csv, show_col_types = FALSE, progress = FALSE)

required_mk <- c(
  "GP2ID", "study", "study_type",
  "wgs", "wgs_prune_reason", "wgs_label",
  "baseline_GP2_phenotype",
  "age_at_sample_collection", "age_of_onset"
)
missing_mk <- setdiff(required_mk, names(mk))
if (length(missing_mk) > 0) {
  stop("Master key missing required columns: ", paste(missing_mk, collapse = ", "))
}

mk <- mk %>%
  mutate(
    study_type = clean_chr(study_type),
    wgs_label  = clean_chr(wgs_label),
    wgs_prune_reason = clean_chr(wgs_prune_reason),
    baseline_GP2_phenotype = clean_chr(baseline_GP2_phenotype),
    study = clean_chr(study)
  )

selected_study_types <- c("Case(/Control)", "Monogenic", "Brain Bank")

mk0 <- mk %>%
  filter(
    wgs == 1,
    is.na(wgs_prune_reason),  # after clean_chr, "" -> NA
    study_type %in% selected_study_types
  )

cat("Total WGS usable samples in selected study types:\n")
print(mk0 %>% count(study_type, name = "n"))

# ------------------------------------------------------------------
# 2) Phenotype recode (baseline_GP2_phenotype only)
# ------------------------------------------------------------------
mk1 <- mk0 %>%
  mutate(
    PHENO_PLINK = case_when(
      baseline_GP2_phenotype %in% c("PD", "Affected_PD") ~ 2L,
      baseline_GP2_phenotype %in% c("Control", "Unaffected") ~ 1L,
      TRUE ~ 0L
    ),
    PHENO = case_when(
      PHENO_PLINK == 2L ~ "Case",
      PHENO_PLINK == 1L ~ "Control",
      TRUE ~ "Other"
    ),
    age_onset_or_recruitment = case_when(
      PHENO_PLINK == 2L ~ coalesce(age_of_onset, age_at_sample_collection),
      PHENO_PLINK == 1L ~ age_at_sample_collection,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(PHENO_PLINK %in% c(1L, 2L))

cat("\nSamples with PD/control phenotype after recoding (before related filter):\n")
print(mk1 %>% count(study_type, PHENO, name = "n") %>% arrange(study_type, PHENO))

# Drop missing ancestry label
na_anc <- mk1 %>% filter(is.na(wgs_label)) %>% nrow()
if (na_anc > 0) cat("\nNOTE: Dropping ", na_anc, " PD/control samples with missing wgs_label.\n", sep = "")
mk1 <- mk1 %>% filter(!is.na(wgs_label))

ancestries <- sort(unique(mk1$wgs_label))
cat("\nAncestries present in WGS PD/control set:", paste(ancestries, collapse = ", "), "\n\n")

# ------------------------------------------------------------------
# 3) Relatedness pruning per ancestry (component-based, using all edges in .related)
# ------------------------------------------------------------------
set.seed(12345)

kept_by_anc <- list()
qc_prune <- list()
qc_sel_rows <- list()

for (anc in ancestries) {
  cat("=== Processing ancestry:", anc, "===\n")
  dat_anc <- mk1 %>% filter(wgs_label == anc)

  cat("  Samples before related filter:", nrow(dat_anc), "\n")

  rel_file <- file.path(related_base, anc, paste0(anc, "_release10.related"))
  res <- infer_unrelated_selection(dat_anc, rel_file)

  keep_ids <- res$keep_ids

  cat("  Related edges used:", res$n_edges, "\n")
  cat("  Clusters (size>1):", res$n_clusters_gt1, "\n")
  cat("  Kept as unrelated :", length(keep_ids), "\n\n")

  kept_by_anc[[anc]] <- keep_ids

  qc_prune[[anc]] <- tibble(
    ancestry = anc,
    n_before = nrow(dat_anc),
    n_edges = res$n_edges,
    n_clusters_gt1 = res$n_clusters_gt1,
    n_kept = length(keep_ids),
    n_dropped = nrow(dat_anc) - length(keep_ids),
    related_file_exists = file.exists(rel_file)
  )

  qc_sel_rows[[anc]] <- res$selection_tbl %>%
    mutate(ancestry = anc)
}

qc_prune_tbl <- bind_rows(qc_prune) %>% arrange(ancestry)
qc_selection <- bind_rows(qc_sel_rows) %>% arrange(ancestry, desc(KEEP_UNRELATED), cluster_id, GP2ID)

write_tsv(qc_prune_tbl, qc_prune_summary)
write_tsv(qc_selection, qc_selection_tbl)

cat("Wrote pruning QC summary to:\n  ", qc_prune_summary, "\n", sep = "")
cat("Wrote pre-PSAM selection table to:\n  ", qc_selection_tbl, "\n\n", sep = "")

mk_unrel_prepsam <- mk1 %>%
  mutate(KEEP_UNRELATED = (GP2ID %in% unlist(kept_by_anc))) %>%
  filter(KEEP_UNRELATED)

cat("Total unrelated samples (before PSAM intersection):", nrow(mk_unrel_prepsam), "\n\n")

# ------------------------------------------------------------------
# 4) Intersect with WGS PSAM (IID-only); write PLINK keep (FID from PSAM)
# ------------------------------------------------------------------
keep_rows <- list()
psam_used <- list()

for (anc in ancestries) {
  psam_obj <- read_psam(anc)
  if (is.null(psam_obj)) {
    warning("No PSAM found for ancestry ", anc, "; skipping PSAM intersection for this ancestry.")
    next
  }

  psam <- psam_obj$data
  psam_path <- psam_obj$path

  ids <- kept_by_anc[[anc]]
  if (length(ids) == 0) next

  psam_sub <- psam %>%
    filter(IID %in% ids) %>%
    select(FID, IID)

  keep_rows[[anc]] <- psam_sub
  psam_used[[anc]] <- tibble(ancestry = anc, psam_path = psam_path, n_matched = nrow(psam_sub))
}

keep_tbl <- bind_rows(keep_rows) %>% distinct()
psam_used_tbl <- bind_rows(psam_used) %>% arrange(ancestry)

write_tsv(psam_used_tbl, qc_psam_used)

cat("Total unrelated samples after PSAM intersection (IID-only match):", nrow(keep_tbl), "\n")
cat("Wrote PSAM-used QC to:\n  ", qc_psam_used, "\n\n", sep = "")

# IMPORTANT: keep file should have NO header for PLINK --keep
write_tsv(keep_tbl, out_keep, col_names = FALSE)
write_lines(sort(unique(keep_tbl$IID)), out_iids)

# ------------------------------------------------------------------
# 5) Final phenotype table aligned to keep_tbl IIDs
# ------------------------------------------------------------------
mk_final <- mk_unrel_prepsam %>%
  inner_join(keep_tbl, by = c("GP2ID" = "IID")) %>%
  mutate(FID = as.character(FID)) %>%
  relocate(any_of(c(
    "study", "FID", "GP2ID",
    "wgs_label", "study_type",
    "baseline_GP2_phenotype",
    "PHENO", "PHENO_PLINK",
    "age_onset_or_recruitment"
  )))

write_tsv(mk_final, out_pheno)

# ------------------------------------------------------------------
# 6) Counts by ancestry
# ------------------------------------------------------------------
counts <- mk_final %>%
  group_by(wgs_label) %>%
  summarise(
    ancestry   = first(wgs_label),
    n_total    = n(),
    n_cases    = sum(PHENO_PLINK == 2L),
    n_controls = sum(PHENO_PLINK == 1L),
    .groups = "drop"
  ) %>%
  select(ancestry, n_total, n_cases, n_controls) %>%
  arrange(ancestry)

counts_all <- counts %>%
  summarise(
    ancestry = "ALL",
    n_total = sum(n_total),
    n_cases = sum(n_cases),
    n_controls = sum(n_controls)
  )

counts_out <- bind_rows(counts, counts_all)
write_tsv(counts_out, out_counts)

cat("\nFinal unrelated case/control counts by ancestry (after PSAM; IID-only match):\n")
print(counts_out)

cat("\nWrote phenotype to:\n  ", out_pheno, "\n", sep = "")
cat("Wrote plink keep file to:\n  ", out_keep, "\n", sep = "")
cat("Wrote IID-only list to:\n  ", out_iids, "\n", sep = "")
cat("Wrote ancestry counts to:\n  ", out_counts, "\n", sep = "")

