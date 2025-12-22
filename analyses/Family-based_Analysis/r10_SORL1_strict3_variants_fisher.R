#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
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

# -----------------------------
# Paths (script-local defaults; env overrides allowed)
# -----------------------------
cc_dir <- Sys.getenv("CC_DIR", unset = file.path(script_dir, "cc_data"))
geno_dir <- Sys.getenv("GENO_DIR", unset = file.path(cc_dir, "plink_SORL1_strict3_variants"))

pheno_path <- Sys.getenv(
  "PHENO_PATH",
  unset = file.path(cc_dir, "r10_wgs_cc_mono_plus_brainbank_UNRELATED_phenotype.tsv")
)

raw_paths <- list(
  EAS = file.path(geno_dir, "chr11_EAS_SORL1_strict3_unrelated.raw"),
  EUR = file.path(geno_dir, "chr11_EUR_SORL1_strict3_unrelated.raw")
)

out_path <- Sys.getenv(
  "OUT_PATH",
  unset = file.path(cc_dir, "r10_SORL1_strict3_variants_carrier_counts_fisher.tsv")
)

dir.create(cc_dir, showWarnings = FALSE, recursive = TRUE)

cat("CC_DIR     :", cc_dir, "\n")
cat("GENO_DIR   :", geno_dir, "\n")
cat("PHENO_PATH :", pheno_path, "\n")
cat("OUT_PATH   :", out_path, "\n\n")

# -----------------------------
# Helpers
# -----------------------------
fmt3 <- function(x) {
  ifelse(
    is.na(x), NA_character_,
    ifelse(
      x == 0, "0",
      ifelse(abs(x) < 0.001,
             format(x, scientific = TRUE, digits = 3),
             format(round(x, 3), nsmall = 3, trim = TRUE))
    )
  )
}

parse_variant_id <- function(variant_col) {
  # e.g. "chr11:121514222:A:C_A" -> "chr11:121514222:A:C"
  sub("_[ACGT]+$", "", variant_col)
}

parse_suffix_allele <- function(variant_col) {
  # e.g. "chr11:..._A" -> "A"
  sub(".*_", "", variant_col)
}

parse_ref_alt <- function(variant_id) {
  # "chr11:pos:REF:ALT" -> c(ref, alt)
  parts <- strsplit(variant_id, ":", fixed = TRUE)[[1]]
  if (length(parts) < 4) return(c(NA_character_, NA_character_))
  c(parts[3], parts[4])
}

is_carrier_from_raw <- function(geno_dosage, variant_col, variant_id) {
  # PLINK .raw columns are allele dosages for the allele in the suffix (e.g. _A)
  # If suffix allele == REF: carrier if dosage < 2 (het=1 or hom-alt=0)
  # If suffix allele == ALT: carrier if dosage > 0
  # Missing may be NA or -9

  if (is.na(geno_dosage)) return(NA)
  if (!is.na(geno_dosage) && geno_dosage == -9) return(NA)

  suffix <- parse_suffix_allele(variant_col)
  ra <- parse_ref_alt(variant_id)
  ref <- ra[1]; alt <- ra[2]

  if (!is.na(ref) && suffix == ref) {
    return(geno_dosage < 2)
  } else if (!is.na(alt) && suffix == alt) {
    return(geno_dosage > 0)
  } else {
    # Fallback: assume suffix is REF-like dosage
    return(geno_dosage < 2)
  }
}

# -----------------------------
# Preconditions
# -----------------------------
if (!file.exists(pheno_path)) stop("Phenotype file missing: ", pheno_path)
for (nm in names(raw_paths)) {
  if (!file.exists(raw_paths[[nm]])) stop("Missing .raw for ", nm, ": ", raw_paths[[nm]])
}

# -----------------------------
# Load phenotype table (already unrelated + PSAM-intersected)
# -----------------------------
pheno <- read_tsv(pheno_path, col_types = cols(), progress = FALSE)

req_cols <- c("GP2ID", "PHENO_PLINK", "wgs_label")
missing <- setdiff(req_cols, names(pheno))
if (length(missing) > 0) stop("Phenotype table missing required cols: ", paste(missing, collapse = ", "))

pheno <- pheno %>%
  transmute(
    IID = as.character(GP2ID),
    ancestry = as.character(wgs_label),
    PHENO_PLINK = suppressWarnings(as.integer(PHENO_PLINK))
  ) %>%
  filter(PHENO_PLINK %in% c(1L, 2L))

cat("Total unrelated CC+Monogenic+BrainBank with PD/control phenotype:", nrow(pheno), "\n\n")

# -----------------------------
# Process each ancestry raw file
# -----------------------------
all_res <- list()

for (anc in names(raw_paths)) {
  raw_path <- raw_paths[[anc]]

  cat("=== Processing ancestry:", anc, "===\n")
  raw <- read_table(raw_path, col_types = cols(), progress = FALSE)

  if (!("IID" %in% names(raw))) stop("No IID column found in .raw: ", raw_path)

  # Genotype columns are chr11:* (include suffix allele)
  geno_cols <- grep("^chr11:", names(raw), value = TRUE)
  if (length(geno_cols) == 0) stop("No genotype columns found in: ", raw_path)

  cat("  Genotype columns in .raw:", paste(geno_cols, collapse = ", "), "\n")

  ph_anc <- pheno %>% filter(ancestry == anc)
  cat("  Phenotype rows for ancestry:", nrow(ph_anc), "\n")

  dat <- raw %>%
    transmute(
      IID = as.character(IID),
      across(all_of(geno_cols), ~ suppressWarnings(as.numeric(.x)))
    ) %>%
    inner_join(ph_anc, by = "IID")

  cat("  Samples after merge with phenotype:", nrow(dat), "\n")

  if (nrow(dat) == 0) {
    warning("No merged samples for ancestry ", anc, " — skipping.")
    cat("\n")
    next
  }

  for (vc in geno_cols) {
    vid <- parse_variant_id(vc)

    carrier_vec <- vapply(
      dat[[vc]],
      FUN = function(g) is_carrier_from_raw(g, vc, vid),
      FUN.VALUE = logical(1)
    )

    nonmissing <- !is.na(carrier_vec)

    n_case <- sum(dat$PHENO_PLINK == 2L & nonmissing)
    n_ctrl <- sum(dat$PHENO_PLINK == 1L & nonmissing)

    case_carriers <- sum(dat$PHENO_PLINK == 2L & carrier_vec, na.rm = TRUE)
    ctrl_carriers <- sum(dat$PHENO_PLINK == 1L & carrier_vec, na.rm = TRUE)

    case_noncarriers <- n_case - case_carriers
    ctrl_noncarriers <- n_ctrl - ctrl_carriers

    mat <- matrix(
      c(case_carriers, case_noncarriers,
        ctrl_carriers, ctrl_noncarriers),
      nrow = 2, byrow = TRUE
    )

    ft <- suppressWarnings(fisher.test(mat, alternative = "greater"))

    or_est <- unname(ft$estimate)
    if (length(or_est) == 0) or_est <- NA_real_

    case_freq <- ifelse(n_case > 0, case_carriers / n_case, NA_real_)
    ctrl_freq <- ifelse(n_ctrl > 0, ctrl_carriers / n_ctrl, NA_real_)

    all_res[[length(all_res) + 1]] <- tibble(
      ancestry = anc,
      variant_col = vc,
      variant_id = vid,

      n_case = n_case,
      n_ctrl = n_ctrl,

      case_carriers = case_carriers,
      case_noncarriers = case_noncarriers,
      ctrl_carriers = ctrl_carriers,
      ctrl_noncarriers = ctrl_noncarriers,

      case_freq = case_freq,
      ctrl_freq = ctrl_freq,

      case_carriers_over_total = paste0(case_carriers, "/", n_case),
      ctrl_carriers_over_total = paste0(ctrl_carriers, "/", n_ctrl),

      case_freq_fmt = fmt3(case_freq),
      ctrl_freq_fmt = fmt3(ctrl_freq),

      fisher_p_one_sided = ft$p.value,
      fisher_OR = as.numeric(or_est),

      fisher_p_fmt = ifelse(
        is.na(ft$p.value), NA_character_,
        ifelse(ft$p.value < 0.001,
               format(ft$p.value, scientific = TRUE, digits = 3),
               format(round(ft$p.value, 3), nsmall = 3, trim = TRUE))
      )
    )
  }

  cat("  Completed ancestry:", anc, "\n\n")
}

res <- bind_rows(all_res) %>%
  arrange(ancestry, variant_id)

write_tsv(res, out_path)

cat("Wrote carrier counts + Fisher tests to:\n  ", out_path, "\n\n", sep = "")
print(res %>% select(
  ancestry, variant_id, n_case, n_ctrl,
  case_carriers, ctrl_carriers,
  case_freq_fmt, ctrl_freq_fmt,
  fisher_p_fmt, fisher_OR
))

