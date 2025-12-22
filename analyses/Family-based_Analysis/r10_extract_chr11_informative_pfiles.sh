#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Purpose:
#   Extract chr11 SORL1 region pfiles from GP2 r10 joint-called WGS
#   PLINK2 pfiles, restricted to "informative WGS" individuals per ancestry.
#
# Output location:
#   By default, writes outputs next to this script (script-local).
#   You can override with OUT_DIR / OUT_BASE_DIR if needed.
#
# Env overrides:
#   PLINK2, WGS_PLINK_BASE, OUT_DIR, OUT_BASE_DIR
#   SORL1_CHR, SORL1_START, SORL1_END, THREADS
#   ANCESTRIES (space-separated list) to override auto-detection
# ------------------------------------------------------------

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Where keep files live (default: next to this script)
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}}"

PLINK2="${PLINK2:-/home/jupyter/bernabe/programs/plink2}"
WGS_PLINK_BASE="${WGS_PLINK_BASE:-/home/jupyter/workspace/gp2_tier2_eu_release10/wgs/deepvariant_joint_calling/plink}"

OUT_BASE_DIR="${OUT_BASE_DIR:-${OUT_DIR}/plink_chr11_informative}"
mkdir -p "${OUT_BASE_DIR}"

SORL1_CHR="${SORL1_CHR:-11}"
SORL1_START="${SORL1_START:-121452314}"
SORL1_END="${SORL1_END:-121633763}"

THREADS="${THREADS:-4}"

# Sanity checks
if [[ ! -x "${PLINK2}" ]]; then
  echo "ERROR: PLINK2 not found or not executable at: ${PLINK2}" >&2
  exit 1
fi
if [[ ! -d "${WGS_PLINK_BASE}" ]]; then
  echo "ERROR: WGS_PLINK_BASE directory not found: ${WGS_PLINK_BASE}" >&2
  exit 1
fi
if [[ ! -d "${OUT_DIR}" ]]; then
  echo "ERROR: OUT_DIR not found: ${OUT_DIR}" >&2
  exit 1
fi

echo "SCRIPT_DIR    : ${SCRIPT_DIR}"
echo "OUT_DIR       : ${OUT_DIR}"
echo "OUT_BASE_DIR  : ${OUT_BASE_DIR}"
echo "PLINK2        : ${PLINK2}"
echo "WGS_PLINK_BASE: ${WGS_PLINK_BASE}"
echo "Region        : chr${SORL1_CHR}:${SORL1_START}-${SORL1_END}"
echo "THREADS       : ${THREADS}"
echo

# Determine ancestries:
# - If ANCESTRIES env var is set, use it (space-separated).
# - Otherwise auto-detect from *keep.tsv files in OUT_DIR.
declare -a ANC_ARR=()
if [[ -n "${ANCESTRIES:-}" ]]; then
  read -r -a ANC_ARR <<< "${ANCESTRIES}"
else
  mapfile -t ANC_ARR < <(
    cd "${OUT_DIR}" && \
      ls *_r10_informative_keep.tsv 2>/dev/null | \
      sed 's/_r10_informative_keep\.tsv$//' | \
      grep -v '^ALL$' | \
      sort -u
  )
fi

if [[ "${#ANC_ARR[@]}" -eq 0 ]]; then
  echo "ERROR: No ancestries detected. Expected *_r10_informative_keep.tsv in: ${OUT_DIR}" >&2
  exit 1
fi

echo "ANCESTRIES to process: ${ANC_ARR[*]}"
echo

for ANC in "${ANC_ARR[@]}"; do
  KEEP_FILE="${OUT_DIR}/${ANC}_r10_informative_keep.tsv"
  if [[ ! -f "${KEEP_FILE}" ]]; then
    echo "[$ANC] Keep file not found at ${KEEP_FILE} – skipping."
    continue
  fi

  IN_BASE="${WGS_PLINK_BASE}/${ANC}/chr11_${ANC}_release10"
  if [[ ! -f "${IN_BASE}.pgen" || ! -f "${IN_BASE}.pvar" || ! -f "${IN_BASE}.psam" ]]; then
    echo "[$ANC] chr11 pfile not found at ${IN_BASE}.{pgen,pvar,psam} – skipping."
    continue
  fi

  ANC_OUT_DIR="${OUT_BASE_DIR}/${ANC}"
  mkdir -p "${ANC_OUT_DIR}"

  OUT_BASE="${ANC_OUT_DIR}/chr11_${ANC}_r10_informative_SORL1"

  echo "[$ANC] Extracting informative chr11 SORL1 region -> ${OUT_BASE}.{pgen,pvar,psam}"
  "${PLINK2}" \
    --pfile "${IN_BASE}" \
    --chr "${SORL1_CHR}" \
    --from-bp "${SORL1_START}" \
    --to-bp "${SORL1_END}" \
    --keep "${KEEP_FILE}" \
    --threads "${THREADS}" \
    --make-pgen \
    --out "${OUT_BASE}"

  echo "[$ANC] Done."
  echo
done

echo "All requested ancestries processed."

