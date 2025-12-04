* **Project:** ADRD-SORL1-Biobanks
* **Version:**R/4.3.1
* **Last Updated:** 14-Jun-2025

## Notebook Overview
#Fine mapping

# Install required packages 
install.packages("data.table")
install.packages("robustbase")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("devtools")
install.packages("dplyr")

# Load libraries
library(data.table)
library(robustbase)
library(tidyr)
library(ggplot2)
library(devtools)
library(dplyr)

# Install and load coloc package
devtools::install_github("chr1swallace/coloc")
library(coloc)

# Load dataset
dataset1 <- fread("~/Desktop/CH_GWAS_formatted_input.tsv", header = TRUE, sep = "\t")

# Add required metadata
dataset1$type <- "cc"               # case-control study
dataset1$s <-   0.4857             # proportion of cases
dataset1$N <-  2240               # total sample size

# Compute squared standard error
dataset_final <- dataset1 %>%
  mutate(StdErr_squared = StdErr^2)

# Prepare output with required columns and new names
output <- dataset_final[, .(
  SNP = MarkerName,
  beta = Effect,
  P = `P-value`,
  varbeta = StdErr_squared,
  type,
  s,
  N
)]

# Drop rows with missing SNPs and duplications
output <- output[!is.na(SNP)]
output <- output[!is.na(beta)]
output <- output[!duplicated(output$SNP), ]

# Save formatted input to Desktop
fwrite(output, file = "~/Desktop/CH_GWAS.csv", na = "NA", quote = FALSE, row.names = FALSE, sep = "\t")

# Run fine-mapping
results <- finemap.abf(dataset = list(
  snp = output$SNP,
  beta = output$beta,
  varbeta = output$varbeta,
  N = output$N[1],
  s = output$s[1],
  type = output$type[1]
))

# Check row counts
cat("Rows in results:", nrow(results), "\n")
cat("Rows in output:", nrow(output), "\n")

# Ensure that results and output have the same number of rows
if (nrow(results) == nrow(output)) {
  # Combine results with input data
  combo <- cbind(results, output)
} else {
  # If rows differ, merge based on SNPs to align the rows
  combo <- merge(results, output, by.x = "snp", by.y = "SNP", all.x = TRUE)
}

# Filter hits with posterior probability > 0.6
hits <- subset(combo, SNP.PP > 0.6)

# Save full results to Desktop
fwrite(combo, file = "~/Desktop/CH_GWAS_fine_map.csv", na = "NA", quote = FALSE, row.names = FALSE, sep = ",")

# Optionally: save filtered hits with posterior probability > 0.6 to Desktop
fwrite(hits, file = "~/Desktop/CH_GWAS_hits_PP_gt_0.6.csv", na = "NA", quote = FALSE, row.names = FALSE, sep = ",")
