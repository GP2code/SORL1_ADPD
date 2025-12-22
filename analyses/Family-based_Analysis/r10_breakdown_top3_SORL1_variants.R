#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

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

# -----------------------------
# Inputs (script-local defaults)
# -----------------------------
seg_fam_path <- Sys.getenv("SEG_FAM_PATH", unset = file.path(out_dir, "r10_SORL1_segregation_by_family.tsv"))
ped_path     <- Sys.getenv("PED_PATH", unset = file.path(out_dir, "r10_master_ped_informative_families.tsv"))
ranked_path  <- Sys.getenv("RANKED_PATH", unset = file.path(out_dir, "r10_SORL1_variants_all_PD_carriers_ranked.tsv"))

# KING/relatedness sources
related_base <- Sys.getenv(
  "RELATED_BASE",
  unset = "/home/jupyter/workspace/gp2_tier2_eu_release10/wgs/deepvariant_joint_calling/related_samples"
)

# Optional PLINK2 fallback (only used if no pairs found in .related)
run_plink_fallback <- tolower(Sys.getenv("RUN_PLINK_FALLBACK", unset = "true")) %in% c("1","true","t","yes","y")
plink2 <- Sys.getenv("PLINK2", unset = "/home/jupyter/bernabe/programs/plink2")
wgs_plink_base <- Sys.getenv(
  "WGS_PLINK_BASE",
  unset = "/home/jupyter/workspace/gp2_tier2_eu_release10/wgs/deepvariant_joint_calling/plink"
)
fallback_chr <- suppressWarnings(as.integer(Sys.getenv("FALLBACK_CHR", unset = "1")))
if (is.na(fallback_chr) || fallback_chr < 1) fallback_chr <- 1
threads <- suppressWarnings(as.integer(Sys.getenv("THREADS", unset = "16")))
if (is.na(threads) || threads < 1) threads <- 16

top_n <- suppressWarnings(as.integer(Sys.getenv("TOP_N", unset = "3")))
if (is.na(top_n) || top_n < 1) top_n <- 3

# -----------------------------
# Outputs
# -----------------------------
out_breakdown_tsv <- file.path(out_dir, "r10_SORL1_top3_family_breakdown.tsv")
out_rel_tsv       <- file.path(out_dir, "r10_SORL1_top3_family_relationship_summary.tsv")
out_rel_md        <- file.path(out_dir, "r10_SORL1_top3_family_relationship_summary.md")

cat("OUT_DIR           :", out_dir, "\n")
cat("SEG_FAM_PATH      :", seg_fam_path, "\n")
cat("PED_PATH          :", ped_path, "\n")
cat("RANKED_PATH       :", ranked_path, "\n")
cat("RELATED_BASE      :", related_base, "\n")
cat("RUN_PLINK_FALLBACK:", run_plink_fallback, "\n")
cat("PLINK2            :", plink2, "\n")
cat("WGS_PLINK_BASE    :", wgs_plink_base, "\n")
cat("FALLBACK_CHR      :", fallback_chr, "\n")
cat("THREADS           :", threads, "\n")
cat("TOP_N             :", top_n, "\n\n")

stopifnot(file.exists(seg_fam_path), file.exists(ped_path), file.exists(ranked_path))

# -----------------------------
# 1) Load main tables
# -----------------------------
seg_fam <- read_tsv(seg_fam_path, col_types = cols())
ped <- read_tsv(ped_path, col_types = cols())
ranked <- read_tsv(ranked_path, col_types = cols())

req_seg <- c("ancestry","variant_id","family_id","n_nonmissing","n_PD","n_ctrl","n_other",
             "n_carrier_PD","n_carrier_ctrl","n_carrier_other","n_noncarrier_PD")
if (!all(req_seg %in% names(seg_fam))) stop("Seg table missing cols: ", paste(setdiff(req_seg, names(seg_fam)), collapse=", "))

req_ped <- c("family_id","id","parental_id","maternal_id","sex","phenotype","ancestry")
if (!all(req_ped %in% names(ped))) stop("PED missing cols: ", paste(setdiff(req_ped, names(ped)), collapse=", "))

req_rank <- c("ancestry","variant_id")
if (!all(req_rank %in% names(ranked))) stop("Ranked missing cols: ", paste(setdiff(req_rank, names(ranked)), collapse=", "))

seg_fam <- seg_fam %>%
  mutate(
    ancestry = as.character(ancestry),
    variant_id = as.character(variant_id),
    family_id = as.character(family_id),
    across(c(n_nonmissing,n_PD,n_ctrl,n_other,n_carrier_PD,n_carrier_ctrl,n_carrier_other,n_noncarrier_PD),
           ~ suppressWarnings(as.integer(.x)))
  )

ped <- ped %>%
  mutate(
    family_id = as.character(family_id),
    id = as.character(id),
    parental_id = as.character(parental_id),
    maternal_id = as.character(maternal_id),
    ancestry = as.character(ancestry),
    phenotype_i = suppressWarnings(as.integer(phenotype))
  )

ranked_top <- ranked %>%
  mutate(ancestry = as.character(ancestry), variant_id = as.character(variant_id)) %>%
  slice_head(n = top_n) %>%
  select(ancestry, variant_id) %>%
  distinct()

cat("Selected variants (TOP_N from ranked):\n")
print(ranked_top)
cat("\n")

# -----------------------------
# 2) Restrict to supporting families only
# -----------------------------
seg_sel <- seg_fam %>%
  inner_join(ranked_top, by = c("ancestry","variant_id")) %>%
  mutate(fam_all_PD_carriers = (n_PD > 0 & n_carrier_PD > 0 & n_noncarrier_PD == 0)) %>%
  filter(fam_all_PD_carriers)

if (nrow(seg_sel) == 0) stop("No family-variant rows meet the all-PD-carriers criterion. Nothing to do.")

fam_list <- seg_sel %>% distinct(family_id, ancestry)

cat("Families to annotate:", nrow(fam_list), "\n\n")

# -----------------------------
# 3) PED-based hints
# -----------------------------
ped_family_hints <- function(df) {
  ids <- df$id
  has_parent_child_edge <- any(df$parental_id %in% ids | df$maternal_id %in% ids)

  shared_moms <- df %>% filter(!is.na(maternal_id), maternal_id != "0") %>% count(maternal_id) %>% filter(n >= 2)
  shared_dads <- df %>% filter(!is.na(parental_id), parental_id != "0") %>% count(parental_id) %>% filter(n >= 2)

  list(
    n_members = nrow(df),
    has_parent_child_edge = has_parent_child_edge,
    has_shared_mother = nrow(shared_moms) > 0,
    has_shared_father = nrow(shared_dads) > 0,
    shared_mother_id = if (nrow(shared_moms) > 0) shared_moms$maternal_id[1] else NA_character_,
    shared_father_id = if (nrow(shared_dads) > 0) shared_dads$parental_id[1] else NA_character_
  )
}

# -----------------------------
# 4) Load ancestry .related (KING) and extract within-family pairs
# -----------------------------
read_related_file <- function(anc) {
  rel_file <- file.path(related_base, anc, paste0(anc, "_release10.related"))
  if (!file.exists(rel_file)) return(NULL)

  # No header; comment lines begin with '#'
  read_delim(
    rel_file,
    delim = ",",
    comment = "#",
    col_names = c("FID1","IID1","FID2","IID2","NSNP","HETHET","IBS0","KINSHIP","REL"),
    col_types = cols(
      FID1 = col_character(), IID1 = col_character(),
      FID2 = col_character(), IID2 = col_character(),
      NSNP = col_double(), HETHET = col_double(),
      IBS0 = col_double(), KINSHIP = col_double(),
      REL = col_character()
    )
  )
}

# Cache related tables by ancestry
rel_cache <- new.env(parent = emptyenv())

get_family_pairs_from_related <- function(anc, ids) {
  if (!exists(anc, envir = rel_cache, inherits = FALSE)) {
    assign(anc, read_related_file(anc), envir = rel_cache)
  }
  rel <- get(anc, envir = rel_cache, inherits = FALSE)
  if (is.null(rel)) return(tibble())

  rel %>%
    filter(IID1 %in% ids, IID2 %in% ids) %>%
    transmute(
      source = "related_file",
      IID1, IID2,
      IBS0 = as.numeric(IBS0),
      KINSHIP = as.numeric(KINSHIP)
    ) %>%
    distinct()
}

# -----------------------------
# 5) Optional PLINK2 KING fallback (writes .kin0 and reads it)
# -----------------------------
read_kin0 <- function(path) {
  if (!file.exists(path)) return(tibble())
  x <- read_table(path, col_types = cols())
  if ("#FID1" %in% names(x)) x <- x %>% rename(FID1 = `#FID1`)
  x %>%
    transmute(
      source = "plink2_king_chr",
      IID1 = as.character(IID1),
      IID2 = as.character(IID2),
      IBS0 = as.numeric(IBS0),
      KINSHIP = as.numeric(KINSHIP)
    )
}

run_plink2_king_for_family <- function(family_id, anc, ids) {
  if (!run_plink_fallback) return(tibble())

  if (!file.exists(plink2)) {
    cat("  ! PLINK2 not found at", plink2, "- cannot run fallback for", family_id, "\n")
    return(tibble())
  }

  pfile <- file.path(wgs_plink_base, anc, paste0("chr", fallback_chr, "_", anc, "_release10"))
  if (!file.exists(paste0(pfile, ".pgen"))) {
    cat("  ! PFILE not found for fallback:", pfile, "\n")
    return(tibble())
  }

  keep_path <- file.path(out_dir, paste0("qc_", family_id, "_", anc, "_chr", fallback_chr, "_king.keep"))
  out_prefix <- file.path(out_dir, paste0("qc_", family_id, "_", anc, "_chr", fallback_chr, "_king"))
  kin0_path <- paste0(out_prefix, ".kin0")

  if (file.exists(kin0_path)) {
    return(read_kin0(kin0_path))
  }

  # Write keep file (FID IID)
  keep_df <- tibble(FID = family_id, IID = ids)
  write_tsv(keep_df, keep_path, col_names = FALSE)

  cat("  * Running PLINK2 KING fallback for", family_id, "(", anc, ") on chr", fallback_chr, "...\n", sep = "")
  args <- c(
    "--pfile", pfile,
    "--keep", keep_path,
    "--make-king-table", "rel-check",
    "--threads", as.character(threads),
    "--out", out_prefix
  )

  rc <- tryCatch({
    system2(plink2, args = args)
  }, error = function(e) 1)

  if (!file.exists(kin0_path)) {
    cat("  ! Fallback did not produce", kin0_path, "for", family_id, "\n")
    return(tibble())
  }

  read_kin0(kin0_path)
}

# -----------------------------
# 6) Pairwise -> relationship calls
# -----------------------------
ibs0_po_thresh <- 1e-4

classify_pair <- function(kinship, ibs0) {
  if (is.na(kinship) || is.na(ibs0)) return(NA_character_)
  if (kinship >= 0.20) {
    if (ibs0 <= ibs0_po_thresh) return("parent_offspring_like")
    return("full_sib_like")
  }
  if (kinship >= 0.088) return("second_degree_like")
  if (kinship >= 0.044) return("third_degree_like")
  "unrelated_or_distant"
}

# -----------------------------
# 7) Build family-level relationship summary
# -----------------------------
family_relationship_summary <- fam_list %>%
  pmap_dfr(function(family_id, ancestry) {
    ped_fam <- ped %>% filter(family_id == !!family_id, ancestry == !!ancestry)
    if (nrow(ped_fam) == 0) {
      return(tibble(
        family_id = family_id, ancestry = ancestry,
        n_members = NA_integer_,
        relationship_class = "missing_ped",
        relationship_summary = "missing PED rows for family/ancestry",
        evidence = NA_character_
      ))
    }

    hints <- ped_family_hints(ped_fam)
    ids <- ped_fam$id

    # Try .related first
    pairs <- get_family_pairs_from_related(ancestry, ids)

    # Fallback if no pairs found
    if (nrow(pairs) == 0 && length(ids) >= 2) {
      pairs <- run_plink2_king_for_family(family_id, ancestry, ids)
    }

    # Pair-level classifications
    if (nrow(pairs) > 0) {
      pairs <- pairs %>%
        mutate(pair_class = mapply(classify_pair, KINSHIP, IBS0))
      kin_mean <- mean(pairs$KINSHIP, na.rm = TRUE)
      kin_max  <- max(pairs$KINSHIP, na.rm = TRUE)
      ibs0_min <- min(pairs$IBS0, na.rm = TRUE)
      ibs0_max <- max(pairs$IBS0, na.rm = TRUE)

      n_po <- sum(pairs$pair_class == "parent_offspring_like", na.rm = TRUE)
      n_fs <- sum(pairs$pair_class == "full_sib_like", na.rm = TRUE)
      n_2d <- sum(pairs$pair_class == "second_degree_like", na.rm = TRUE)

      rel_source <- paste(unique(pairs$source), collapse = ",")
    } else {
      kin_mean <- NA_real_; kin_max <- NA_real_
      ibs0_min <- NA_real_; ibs0_max <- NA_real_
      n_po <- 0; n_fs <- 0; n_2d <- 0
      rel_source <- "none"
    }

    # Final family-level interpretation
    relationship_class <- "unclear"
    relationship_summary <- "unclear"
    evidence_bits <- c(
      paste0("PED:parent_child_edge=", hints$has_parent_child_edge),
      paste0("PED:shared_mother=", hints$has_shared_mother),
      paste0("PED:shared_father=", hints$has_shared_father),
      paste0("KINGsrc=", rel_source),
      paste0("KIN_mean=", ifelse(is.na(kin_mean), "NA", sprintf("%.6f", kin_mean))),
      paste0("KIN_max=",  ifelse(is.na(kin_max),  "NA", sprintf("%.6f", kin_max))),
      paste0("IBS0_min=", ifelse(is.na(ibs0_min), "NA", sprintf("%.8f", ibs0_min)))
    )

    if (hints$n_members == 2 && hints$has_parent_child_edge) {
      relationship_class <- "parent_child"
      relationship_summary <- "parent–child (PED encodes it; IBS0≈0 supports PO)"
    } else if (hints$n_members == 2 && !is.na(kin_max) && kin_max >= 0.20 && !is.na(ibs0_min) && ibs0_min <= ibs0_po_thresh) {
      relationship_class <- "parent_child_like"
      relationship_summary <- "parent–child-like (IBS0≈0; KINSHIP≈0.25)"
    } else if (hints$n_members >= 2 && hints$has_shared_mother && !is.na(kin_mean) && kin_mean >= 0.20 && n_po == 0 && n_2d == 0) {
      # Shared mother + all pairs look like FS (no PO-like pairs)
      if (hints$n_members == 2) {
        relationship_class <- "siblings"
        relationship_summary <- "siblings (shared mother in PED; IBS0>0 supports sib vs PO)"
      } else {
        relationship_class <- "sibship"
        relationship_summary <- paste0(hints$n_members, " siblings (shared mother in PED; all pairs first-degree)")
      }
    } else if (hints$n_members == 2 && !is.na(kin_mean) && kin_mean >= 0.088 && kin_mean < 0.20) {
      relationship_class <- "second_degree_like"
      relationship_summary <- paste0("second-degree-like (KINSHIP=", sprintf("%.6f", kin_mean),
                                    "; IBS0>0; half-sib/avuncular/grandparent)")
    } else if (hints$n_members == 2 && !is.na(kin_mean) && kin_mean >= 0.20 && !is.na(ibs0_min) && ibs0_min > ibs0_po_thresh) {
      relationship_class <- "siblings_like"
      relationship_summary <- "siblings-like (KINSHIP≈0.25; IBS0>0)"
    } else if (hints$n_members >= 3 && n_po >= 2) {
      relationship_class <- "parent_plus_children_like"
      relationship_summary <- "likely parent + children (multiple PO-like pairs)"
    } else if (hints$n_members >= 2 && hints$has_shared_mother) {
      relationship_class <- "sibship_shared_mother_unconfirmed"
      relationship_summary <- "sibship suggested by PED (shared mother), but KING evidence is incomplete"
    } else if (hints$n_members == 2 && (is.na(kin_mean) || is.na(ibs0_min))) {
      relationship_class <- "unknown_no_king"
      relationship_summary <- "relationship unknown (no .related/king data found)"
    }

    tibble(
      family_id = family_id,
      ancestry = ancestry,
      n_members = hints$n_members,
      relationship_class = relationship_class,
      relationship_summary = relationship_summary,
      king_source = rel_source,
      king_n_pairs = ifelse(exists("pairs") && is.data.frame(pairs), nrow(pairs), 0L),
      king_kinship_mean = kin_mean,
      king_kinship_max = kin_max,
      king_ibs0_min = ibs0_min,
      king_ibs0_max = ibs0_max,
      evidence = paste(evidence_bits, collapse = "; ")
    )
  })

write_tsv(family_relationship_summary, out_rel_tsv)

# Markdown summary (bullet list)
md_lines <- c(
  "# SORL1 top-variant family relationship summary",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC")),
  ""
)
for (i in seq_len(nrow(family_relationship_summary))) {
  r <- family_relationship_summary[i, ]
  md_lines <- c(md_lines, paste0("- **", r$family_id, "**: ", r$relationship_summary))
}
writeLines(md_lines, out_rel_md)

cat("Wrote family relationship summary TSV:\n  ", out_rel_tsv, "\n", sep = "")
cat("Wrote family relationship summary MD:\n  ", out_rel_md, "\n\n", sep = "")

# -----------------------------
# 8) Build enriched breakdown table (variant x family + family structure + relationship summary)
# -----------------------------
fam_struct <- ped %>%
  group_by(family_id, ancestry) %>%
  summarise(
    family_n_members = n(),
    family_n_PD      = sum(phenotype_i == 2, na.rm = TRUE),
    family_n_ctrl    = sum(phenotype_i == 1, na.rm = TRUE),
    family_n_other   = family_n_members - family_n_PD - family_n_ctrl,
    n_founders       = sum(parental_id == "0" & maternal_id == "0"),
    n_with_parentinfo= sum(parental_id != "0" | maternal_id != "0"),
    .groups = "drop"
  )

top_breakdown <- seg_sel %>%
  select(
    ancestry, variant_id, family_id,
    n_nonmissing, n_PD, n_ctrl, n_other,
    n_carrier_PD, n_carrier_ctrl, n_carrier_other,
    n_noncarrier_PD
  ) %>%
  left_join(fam_struct, by = c("family_id","ancestry")) %>%
  left_join(family_relationship_summary, by = c("family_id","ancestry")) %>%
  arrange(ancestry, variant_id, family_id)

write_tsv(top_breakdown, out_breakdown_tsv)

cat("Wrote enriched TOP_N family breakdown to:\n  ", out_breakdown_tsv, "\n", sep = "")

