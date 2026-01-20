
# *SORL1* Across Neurodegenerative Diseases: A Multi-Ancestry Biobank-Scale Assessment

`GP2 ❤️ Open Science 😍`

[![DOI](https://zenodo.org/badge/1028658758.svg)](https://doi.org/10.5281/zenodo.17819222)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Last Updated:** January 2026

---

## **Summary**

This repository contains the code, data workflows, and results associated with the manuscript titled:  
**“Is SORL1 a common genetic target across neurodegenerative diseases?: A multi-ancestry biobank-scale assessment.”**

In this study, we performed the largest global genetic analysis of the *SORL1* gene across Alzheimer’s disease (AD), related dementias (RD), and Parkinson’s disease (PD) using whole-genome sequencing (WGS) and imputed data from six large biobank-scale resources. We analyzed 67,749 cases and 111,969 controls from 11 genetically defined ancestries, evaluating rare protein-altering and splicing variants, performing burden and association tests, and expression impact across diseases and populations. In addition, we conducted a family-based analysis of 398 informative families with at least one PD case and at least two members.

Highlights of this work include:

- Identification of 53 potentially disease-causing *SORL1* variants, including 41 novel variants
- Replication of *SORL1* variants across AD, RD, and PD, highlighting pleiotropic effects
- Discovery of ancestry-specific and cross-ancestry associations
- A family-based analysis identified a rare predicted-damaging variant in East Asians and two in European ancestries that show evidence of segregation in PD families


---

## **Data Statement**

- **Whole-genome sequencing (WGS) and imputed genotype data** were obtained from the following cohorts:
  - [All of Us (AoU) v7](https://www.researchallofus.org/)
  - [Alzheimer’s Disease Sequencing Project (ADSP) v4](https://www.niagads.org/adsp)
  - [Accelerating Medicines Partnership – Parkinson’s Disease (AMP PD) v3](https://amp-pd.org/)
  - [100,000 Genomes Project (100KGP) v18](https://www.genomicsengland.co.uk/)
  - [UK Biobank (UKB) v18.1](https://www.ukbiobank.ac.uk/)
  - [Global Parkinson’s Genetics Program (GP2) releases v9 and v10](https://gp2.org/)

- All cohorts underwent standardized quality control and ancestry inference using the [GenoTools](https://github.com/dvitale199/GenoTools) pipeline except for 100KGP.
- Access to each dataset must be requested from the respective data platforms.

---

## **Citation**

If you use this repository or find it helpful for your research, please cite the corresponding manuscript:

> **Is SORL1 a Common Genetic Target Across Neurodegenerative Diseases?: A Multi-Ancestry Biobank-Scale Assessment**  
> *Khani M, Yeboah SN, Cerquera-Cleves C, Kedmi A, Bustos BI, Grant SM, Akerman SC, Akçimen F, Lee PS, Reyes-Pérez P, Lange LM, Leonard H, Koretsky MJ, Makarious MB, Schneider Z, Rothstein JD, Merchant K, Mencacci NE, Krainc D, Cookson MR, Singleton AB, Bandres-Ciga S; Global Parkinson's Genetics Program (GP2), 2026*  
> *(DOI: pending)*

---

## **Repository Orientation**

```bash
.
├── ADSP/
│   └── 00_ADSP.ipynb
├── AllofUs/
│   └── 00_AllofUs.ipynb
├── AMP-PD/
│   ├── 00_AMP_Clinical data.ipynb
│   ├── 00_AMP_PD_DLB.ipynb
│   └── 00_AMP_PD.ipynb
├── GP2/
│   └── 00_GP2.ipynb
├── UKBiobank/
│   └── 00_UKB.ipynb
├── 06_Fine_mapping.R
├── Visualization_Protein_Structure/
│   └── q92673_model.cif
│   └── label_SORL1_PD.pml
│   └── label_SORL1_ADRD.pml
├── Family-based_Analysis/
│   └──r10_annotate_plink2_fams.sh
│   └──r10_extract_SORL1_rare_damaging_pfiles.sh
│   └──r10_extract_chr11_informative_pfiles.sh
│   └──r10_export_SORL1_rare_damaging_raw.sh
│   └──r10_extract_SORL1_strict3_variants_cc.sh
│   └──r10_make_informative_keep_files.R
│   └──r10_rank_SORL1_all_PD_carrier_variants.R
│   └──r10_SORL1_strict3_variants_fisher.R
│   └──r10_merge_core_fams_with_inferred_fams.R
│   └──r10_mark_informative_families.R
│   └──r10_annotate_core_ped_with_master.R
│   └──r10_kinship_clusters_to_ped.R
│   └──r10_summarize_SORL1_segregation_by_ancestry.R
│   └──r10_filter_SORL1_by_ancestry.R
│   └──r10_rebuild_unrelated_cc_mono_brainbank_IIDpsam.R
│   └──r10_breakdown_top3_SORL1_variants.R
└── LICENSE
```
---

## Key Analyses
1. Variant Filtering
    * Filtering based on CADD > 20, MAC ≥ 2, protein-altering or splicing impact, and case-only presence
2. Association and Burden Analysis
    * Case-control comparisons using logistic regression and burden tests (SKAT-O) via RVTESTS
3. Fine Mapping
    * Bayesian colocalization using coloc with diverse AD/PD GWAS datasets to nominate causal variants
4. eQTL Annotation
     * Multivariate brain eQTLs using the mmQTL resource from PsychENCODE, GTEx, and ROSMAP
5. Family-based analysis
     * Family-based analysis of 398 informative families with at least one PD case and at least two members

## Analysis Notebooks / Scripts
| **Scripts**                 | **Description**                                                                 |
| ----------------------------------- | ------------------------------------------------------------------------------- |
| `06_Fine_mapping.R`                 | Colocalization fine-mapping using multiple AD/PD GWAS datasets across ancestries   |
| `ADSP/00_ADSP.ipynb`                | Variant filtering, association, and burden testing in the ADSP cohort           |
| `AllofUs/00_AllofUs.ipynb`          | *SORL1* variant exploration in AD, RD, and PD within the AoU dataset                 |
| `AMP-PD/00_AMP_PD.ipynb`            | PD-focused analysis of *SORL1* variants from AMP PD WGS data                    |
| `AMP-PD/00_AMP_PD_DLB.ipynb`        | Variant analysis in Dementia with Lewy Bodies (DLB) cases from AMP PD           |
| `AMP-PD/00_AMP_Clinical data.ipynb` | Clinical characterization and cognitive profiling of variant carriers in AMP PD |
| `GP2/00_GP2.ipynb`                  | Variant filtering, Case-control association and burden analysis across ancestries using GP2 imputed data   |
| `UKBiobank/00_UKB.ipynb`            | Variant filtering, Gene-based and single-variant association analysis in UKB AD/RD/PD populations  |
|`q92673_model.cif,label_SORL1_PD.pml,label_SORL1_ADRD.pml`                 |Visualization of Protein Structure|
|`r10_annotate_plink2_fams.sh`            | Annotate SORL1 Region Variants|
|`r10_extract_SORL1_rare_damaging_pfiles.sh`            | Extract SORL1 Rare Damaging files|
|`r10_extract_chr11_informative_pfiles.sh`            | Extract chr11 SORL1 region|
|`r10_export_SORL1_rare_damaging_raw.sh`            | Export SORL1 Rare Damaging Raw Genotypes|
|`r10_extract_SORL1_strict3_variants_cc.sh`            | SORL1 Variant Carrier Analysis|
|`r10_make_informative_keep_files.R`            | Build Informative WGS PLINK Keep Files by Ancestry|
|`r10_rank_SORL1_all_PD_carrier_variants.R`            |Rank SORL1 Variants by All-PD-Carrier Families with Strict Segregation Filter |
|`r10_SORL1_strict3_variants_fisher.R`            | Compute Carrier Counts and Fisher Tests for SORL1  Variants|
|`r10_merge_core_fams_with_inferred_fams.R`            |Merge Core and Kinship-Inferred Pedigrees into Master PED Table |
|`r10_mark_informative_families.R`            | Mark Informative Families and Individuals in Master Pedigree|
|`r10_annotate_core_ped_with_master.R`            |Annotate Core Pedigree with Master Key and Extended Clinical Fields |
|`r10_kinship_clusters_to_ped.R`            |Build Kinship-Inferred Families for Non-Core WGS Samples |
|`r10_summarize_SORL1_segregation_by_ancestry.R`            |Summarize SORL1 Rare Variant Segregation by Ancestry |
|`r10_filter_SORL1_by_ancestry.R`            |Filter Rare Damaging SORL1 Variants by Ancestry |
|`r10_rebuild_unrelated_cc_mono_brainbank_IIDpsam.R`            | Build Unrelated WGS Cohort by Ancestry and Study Type|
|`r10_breakdown_top3_SORL1_variants.R`            |Top SORL1 Variants Family Segregation |

## Software
| Software  | Version   | Resource URL                                                               | RRID              | Notes                                  |
| --------- | --------- | -------------------------------------------------------------------------- | ----------------- | -------------------------------------- |
| R         | 4.3.1     | [r-project.org](https://www.r-project.org/)                                | RRID\:SCR\_001905 | Core data analysis and plotting        |
| PLINK     | 1.9 / 2.0 | [cog-genomics.org/plink/](http://www.cog-genomics.org/plink/)              | RRID\:SCR\_001757 | Genotype extraction and frequency calc |
| RVTESTS   | 2.1.0     | [github.com/zhanxw/rvtests](https://github.com/zhanxw/rvtests)             | RRID\:SCR\_019033 | SKAT-O burden testing                  |
| ANNOVAR   | 2019Oct24 | [wannovar.wglab.org](https://wannovar.wglab.org/)                          | RRID\:SCR\_012821 | Variant annotation                                      |
| GenoTools | -         | [github.com/dvitale199/GenoTools](https://github.com/dvitale199/GenoTools) | -                 | Ancestry prediction and QC             |
| coloc     | 6.0.0     | [CRAN - coloc](https://cran.r-project.org/web/packages/coloc/)             | -                 | Fine-mapping & colocalization          |

---

> *For questions, please open an issue or contact the corresponding author.*
