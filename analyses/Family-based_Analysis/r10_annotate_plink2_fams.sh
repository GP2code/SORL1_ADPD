#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Purpose:
#   For each ancestry-specific extracted SORL1-region pfile set:
#     1) export VCF (bgz) from PLINK2
#     2) annotate with ANNOVAR table_annovar.pl (hg38)
#
# Output location:
#   By default, writes outputs next to this script (script-local).
#
# Env overrides:
#   OUT_DIR, PFILE_BASE_DIR, OUT_BASE_DIR
#   PLINK2, ANNOVAR_DIR, HUMANDB, TABLE_ANNOVAR
#   THREADS_PLINK, THREADS_ANNOVAR
#   ANCESTRIES (space-separated) to override auto-detection
# ------------------------------------------------------------

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}}"

PLINK2="${PLINK2:-/home/jupyter/bernabe/programs/plink2}"

ANNOVAR_DIR="${ANNOVAR_DIR:-/home/jupyter/workspace/ws_files/bernabe/programs/annovar}"
HUMANDB="${HUMANDB:-${ANNOVAR_DIR}/humandb}"

# Correct default in your environment:
DEFAULT_TABLE_ANNOVAR="/home/jupyter/bernabe/programs/table_annovar.pl"
TABLE_ANNOVAR="${TABLE_ANNOVAR:-${DEFAULT_TABLE_ANNOVAR}}"
# Fallback if someone copies this elsewhere:
if [[ ! -f "${TABLE_ANNOVAR}" && -f "${ANNOVAR_DIR}/table_annovar.pl" ]]; then
  TABLE_ANNOVAR="${ANNOVAR_DIR}/table_annovar.pl"
fi

PFILE_BASE_DIR="${PFILE_BASE_DIR:-${OUT_DIR}/plink_chr11_informative}"
OUT_BASE_DIR="${OUT_BASE_DIR:-${OUT_DIR}/annovar_r10}"
mkdir -p "${OUT_BASE_DIR}"

THREADS_PLINK="${THREADS_PLINK:-4}"
THREADS_ANNOVAR="${THREADS_ANNOVAR:-20}"

# Sanity checks
if [[ ! -x "${PLINK2}" ]]; then
  echo "ERROR: PLINK2 not found or not executable at: ${PLINK2}" >&2
  exit 1
fi
if [[ ! -d "${ANNOVAR_DIR}" ]]; then
  echo "ERROR: ANNOVAR_DIR not found: ${ANNOVAR_DIR}" >&2
  exit 1
fi
if [[ ! -d "${HUMANDB}" ]]; then
  echo "ERROR: HUMANDB not found: ${HUMANDB}" >&2
  exit 1
fi
if [[ ! -f "${TABLE_ANNOVAR}" ]]; then
  echo "ERROR: table_annovar.pl not found at: ${TABLE_ANNOVAR}" >&2
  exit 1
fi
if [[ ! -d "${PFILE_BASE_DIR}" ]]; then
  echo "ERROR: PFILE_BASE_DIR not found: ${PFILE_BASE_DIR}" >&2
  exit 1
fi

echo "SCRIPT_DIR     : ${SCRIPT_DIR}"
echo "OUT_DIR        : ${OUT_DIR}"
echo "PFILE_BASE_DIR : ${PFILE_BASE_DIR}"
echo "OUT_BASE_DIR   : ${OUT_BASE_DIR}"
echo "PLINK2         : ${PLINK2}"
echo "ANNOVAR_DIR    : ${ANNOVAR_DIR}"
echo "HUMANDB        : ${HUMANDB}"
echo "TABLE_ANNOVAR  : ${TABLE_ANNOVAR}"
echo "THREADS_PLINK  : ${THREADS_PLINK}"
echo "THREADS_ANNOVAR: ${THREADS_ANNOVAR}"
echo

# Determine ancestries to process:
# - If ANCESTRIES env var set, use it.
# - Else auto-detect from directories under PFILE_BASE_DIR.
declare -a ANC_ARR=()
if [[ -n "${ANCESTRIES:-}" ]]; then
  read -r -a ANC_ARR <<< "${ANCESTRIES}"
else
  mapfile -t ANC_ARR < <(
    find "${PFILE_BASE_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort -u
  )
fi

if [[ "${#ANC_ARR[@]}" -eq 0 ]]; then
  echo "ERROR: No ancestry directories found under: ${PFILE_BASE_DIR}" >&2
  exit 1
fi

echo "ANCESTRIES to process: ${ANC_ARR[*]}"
echo

for ANC in "${ANC_ARR[@]}"; do
  echo "==== ${ANC} ===="

  ANC_PFILE_DIR="${PFILE_BASE_DIR}/${ANC}"
  BASE_IN="${ANC_PFILE_DIR}/chr11_${ANC}_r10_informative_SORL1"

  if [[ ! -f "${BASE_IN}.pgen" || ! -f "${BASE_IN}.pvar" || ! -f "${BASE_IN}.psam" ]]; then
    echo "  ! Skipping ${ANC}: no SORL1 pfile set at ${BASE_IN}.{pgen,pvar,psam}"
    echo
    continue
  fi

  ANC_OUT_DIR="${OUT_BASE_DIR}/${ANC}"
  mkdir -p "${ANC_OUT_DIR}"

  VCF_PREFIX="${ANC_OUT_DIR}/chr11_${ANC}_r10_informative_SORL1"
  VCF="${VCF_PREFIX}.vcf.gz"

  # 1) Export VCF from plink2 pfile
  echo "  * Exporting VCF from PLINK2 pfile..."
  set +e
  "${PLINK2}" \
    --pfile "${BASE_IN}" \
    --recode vcf bgz id-paste=iid \
    --threads "${THREADS_PLINK}" \
    --out "${VCF_PREFIX}"
  rc=$?
  set -e

  if [[ $rc -ne 0 ]]; then
    echo "  ! PLINK2 VCF export failed for ${ANC} (exit=${rc}); skipping ANNOVAR for this ancestry."
    echo
    continue
  fi

  if [[ ! -f "${VCF}" ]]; then
    echo "  ! VCF not created for ${ANC} (${VCF}); skipping ANNOVAR."
    echo
    continue
  fi

  # 2) Annotate with ANNOVAR
  echo "  * Annotating with ANNOVAR (hg38: refGene + gnomAD v4.1 exome/genome + dbnsfp47a)..."
  perl "${TABLE_ANNOVAR}" \
    "${VCF}" "${HUMANDB}" \
    -buildver hg38 \
    -out "${ANC_OUT_DIR}/chr11_${ANC}_r10_informative_SORL1.annotated" \
    -remove \
    -protocol refGene,gnomad41_exome,gnomad41_genome,dbnsfp47a \
    -operation g,f,f,f \
    -nastring . \
    -polish \
    -vcfinput \
    -thread "${THREADS_ANNOVAR}" \
    -otherinfo

  echo "  * Done ${ANC}"
  echo
done

echo "All requested ancestries finished (annotation step)."

