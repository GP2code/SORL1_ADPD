#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Purpose:
#   Export PLINK2 additive .raw genotype matrices for the rare/damaging
#   SORL1 variant sets (per ancestry).
#
# Output location:
#   By default, writes outputs next to the input pfiles in the same folder
#   (script-local pipeline layout under scripts_for_GP2/).
#
# Env overrides:
#   OUT_DIR, PLINK2, PFILE_BASE_DIR, THREADS
#   ANCESTRIES (space-separated list) to override auto-detection
# ------------------------------------------------------------

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}}"

PLINK2="${PLINK2:-/home/jupyter/bernabe/programs/plink2}"
THREADS="${THREADS:-4}"

# Input pfiles (rare/damaging SORL1 per ancestry)
PFILE_BASE_DIR="${PFILE_BASE_DIR:-${OUT_DIR}/plink_chr11_SORL1_rare_damaging}"

# Sanity checks
if [[ ! -x "${PLINK2}" ]]; then
  echo "ERROR: PLINK2 not found or not executable at: ${PLINK2}" >&2
  exit 1
fi
if [[ ! -d "${PFILE_BASE_DIR}" ]]; then
  echo "ERROR: PFILE_BASE_DIR not found: ${PFILE_BASE_DIR}" >&2
  exit 1
fi

echo "SCRIPT_DIR     : ${SCRIPT_DIR}"
echo "OUT_DIR        : ${OUT_DIR}"
echo "PFILE_BASE_DIR : ${PFILE_BASE_DIR}"
echo "PLINK2         : ${PLINK2}"
echo "THREADS        : ${THREADS}"
echo

# Determine ancestries to process:
# - If ANCESTRIES env var set, use it.
# - Else auto-detect from subdirectories under PFILE_BASE_DIR.
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

  ANC_DIR="${PFILE_BASE_DIR}/${ANC}"
  IN_BASE="${ANC_DIR}/chr11_${ANC}_r10_SORL1_rare_damaging"

  if [[ ! -f "${IN_BASE}.pgen" || ! -f "${IN_BASE}.pvar" || ! -f "${IN_BASE}.psam" ]]; then
    echo "  ! Skipping ${ANC}: no rare/damaging pfile at ${IN_BASE}.{pgen,pvar,psam}"
    echo
    continue
  fi

  # Quick QC: number of samples/variants in this pfile
  n_samp="$(awk 'NR>1{c++} END{print c+0}' "${IN_BASE}.psam")"
  n_var="$(awk 'NR>1{c++} END{print c+0}' "${IN_BASE}.pvar")"
  echo "  * Inputs: ${n_samp} samples, ${n_var} variants"

  echo "  * Exporting PLINK2 --export A (.raw) for ${ANC} ..."
  set +e
  "${PLINK2}" \
    --pfile "${IN_BASE}" \
    --export A \
    --threads "${THREADS}" \
    --out "${IN_BASE}"
  rc=$?
  set -e

  if [[ $rc -ne 0 ]]; then
    echo "  ! PLINK2 export failed for ${ANC} (exit=${rc}); continuing."
    echo
    continue
  fi

  if [[ ! -f "${IN_BASE}.raw" ]]; then
    echo "  ! Expected output not found: ${IN_BASE}.raw (continuing)."
    echo
    continue
  fi

  echo "  * Done ${ANC}. Wrote:"
  echo "    ${IN_BASE}.raw"
  echo "    ${IN_BASE}.log"
  echo
done

echo "All requested ancestries processed (exported .raw files)."

