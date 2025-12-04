# *SORL1* Across Neurodegenerative Diseases: A Multi-Ancestry Biobank-Scale Assessment

`GP2 ❤️ Open Science 😍`

\[DOI: pending]

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Last Updated:** December 2025

---

## **Summary**

This repository contains the code, data workflows, and results associated with the manuscript titled:  
**“Is SORL1 a common genetic target across neurodegenerative diseases?: A multi-ancestry biobank-scale assessment.”**

In this study, we performed the largest global genetic analysis of the *SORL1* gene across Alzheimer’s disease (AD), related dementias (RD), and Parkinson’s disease (PD) using whole-genome sequencing (WGS) and imputed data from six large biobank-scale resources. We analyzed 67,749 cases and 111,969 controls from 11 genetically defined ancestries, evaluating rare protein-altering and splicing variants, performing burden and association tests, and characterizing haplotype structure and expression impact across diseases and populations.

Highlights of this work include:

- Identification of 53 potentially disease-causing *SORL1* variants, including 41 novel variants
- Replication of *SORL1* variants across AD, RD, and PD, highlighting pleiotropic effects
- Discovery of ancestry-specific and cross-ancestry associations
- Haplotype and eQTL analyses showing disease-specific architecture and regulatory implications

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
> *Khani M, Yeboah SN, Cerquera-Cleves C, Kedmi A, Grant SM, Akerman SC, Akçimen F, Lee PS, Reyes-Pérez P, Lange LM, Leonard H, Koretsky MJ, Makarious MB, Schneider Z, Merchant K, Rothstein JD, Singleton AB, Bandres-Ciga S; Global Parkinson's Genetics Program (GP2), 2025*  
> *(DOI: pending)*

---

## **Repository Orientation**

```bash
.
├── ADSP/
│   └── 00_ADSP.ipynb
│   └── 00_Haplotype analysis_ADSP.ipynb
│   └── 00_Visualization_ADSP.R
├── AllofUs/
│   └── 00_AllofUs.ipynb
├── AMP-PD/
│   ├── 00_AMP_Clinical data.ipynb
│   ├── 00_AMP_PD_DLB.ipynb
│   └── 00_AMP_PD.ipynb
├── GP2/
│   └── 00_GP2.ipynb
│   └── 01_Haplotype analysis_WGS_GP2.ipynb
│   └── 02_Haplotype analysis_IMP_GP2.ipynb
│   └── 01_Visualization_GP2.R
├── UKBiobank/
│   └── 00_UKB.ipynb
├── 06_Fine_mapping.R
├── Visualization of Protein Structure
│   └── q92673_model.cif
│   └── label_SORL1_PD.pml
│   └── label_SORL1_ADRD.pml
└── LICENSE
```
---

## Key Analyses
1. Variant Filtering
    * Filtering based on CADD > 20, MAC ≥ 2, protein-altering or splicing impact, and case-only presence
2. Association and Burden Analysis
    * Case-control comparisons using logistic regression and burden tests (SKAT-O) via RVTESTS
3. Haplotype Analysis
    * Construction of ±100 kb haplotype blocks surrounding SORL1 by disease and ancestry using common variants
4. Fine Mapping
    * Bayesian colocalization using coloc with diverse AD/PD GWAS datasets to nominate causal variants
5. eQTL Annotation
     * Multivariate brain eQTLs using the mmQTL resource from PsychENCODE, GTEx, and ROSMAP

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
|`ADSP/00_Haplotype analysis_ADSP.ipynb`|Haplotype analysis in ADSP |
|`GP2/01_Haplotype analysis_WGS_GP2.ipynb`| Haplotype analysis in GP2 WGS data|
|`GP2/02_Haplotype analysis_IMP_GP2.ipynb`|Haplotype analysis in GP2 Imputed data|
|`ADSP/00_Visualization_ADSP.R`|Haplotype data visualization in ADSP|
|`GP2/01_Visualization_GP2.R`|Haplotype data visualization in GP2|
|`q92673_model.cif,label_SORL1_PD.pml,label_SORL1_ADRD.pml`                 |Visualization of Protein Structure|


## Software
| Software  | Version   | Resource URL                                                               | RRID              | Notes                                  |
| --------- | --------- | -------------------------------------------------------------------------- | ----------------- | -------------------------------------- |
| R         | 4.3.1     | [r-project.org](https://www.r-project.org/)                                | RRID\:SCR\_001905 | Core data analysis and plotting        |
| PLINK     | 1.9 / 2.0 | [cog-genomics.org/plink/](http://www.cog-genomics.org/plink/)              | RRID\:SCR\_001757 | Genotype extraction and frequency calc |
| RVTESTS   | 2.1.0     | [github.com/zhanxw/rvtests](https://github.com/zhanxw/rvtests)             | RRID\:SCR\_019033 | SKAT-O burden testing                  |
| ANNOVAR   | 2019Oct24 | [wannovar.wglab.org](https://wannovar.wglab.org/)                          | RRID\:SCR\_012821 | Variant annotation                     |
| LDlink    | -         | [ldlink.nih.gov](https://ldlink.nih.gov/)                                  | RRID\:SCR\_014978 | LD calculation                         |
| GenoTools | -         | [github.com/dvitale199/GenoTools](https://github.com/dvitale199/GenoTools) | -                 | Ancestry prediction and QC             |
| coloc     | 6.0.0     | [CRAN - coloc](https://cran.r-project.org/web/packages/coloc/)             | -                 | Fine-mapping & colocalization          |

---

> *For questions, please open an issue or contact the corresponding author.*
