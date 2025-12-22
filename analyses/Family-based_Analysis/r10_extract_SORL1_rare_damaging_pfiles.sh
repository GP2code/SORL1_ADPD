#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Purpose:
#   Subset the ancestry-specific informative SORL1-region pfiles down to the
#   rare/damaging variant set (IDs produced by r10_filter_SORL1_by_ancestry.R).
#
# Output location:
#   By default, writes outputs next to this script (script-local).
#
# Env overrides:
#   OUT_DIR, PLINK2
#   PFILE_BASE_DIR, ANNOVAR_BASE_DIR, OUT_BASE_DIR
#   THREADS
#   ANCESTRIES (space-separated list) to override auto-detection
# ------------------------------------------------------------

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}}"

PLINK2="${PLINK2:-/home/jupyter/bernabe/programs/plink2}"
THREADS="${THREADS:-4}"

PFILE_BASE_DIR="${PFILE_BASE_DIR:-${OUT_DIR}/plink_chr11_informative}"
ANNOVAR_BASE_DIR="${ANNOVAR_BASE_DIR:-${OUT_DIR}/annovar_r10}"
OUT_BASE_DIR="${OUT_BASE_DIR:-${OUT_DIR}/plink_chr11_SORL1_rare_damaging}"
mkdir -p "${OUT_BASE_DIR}"

# Sanity checks
if [[ ! -x "${PLINK2}" ]]; then
  echo "ERROR: PLINK2 not found or not executable at: ${PLINK2}" >&2
  exit 1
fi
if [[ ! -d "${PFILE_BASE_DIR}" ]]; then
  echo "ERROR: PFILE_BASE_DIR not found: ${PFILE_BASE_DIR}" >&2
  exit 1
fi
if [[ ! -d "${ANNOVAR_BASE_DIR}" ]]; then
  echo "ERROR: ANNOVAR_BASE_DIR not found: ${ANNOVAR_BASE_DIR}" >&2
  exit 1
fi

echo "SCRIPT_DIR      : ${SCRIPT_DIR}"
echo "OUT_DIR         : ${OUT_DIR}"
echo "PFILE_BASE_DIR  : ${PFILE_BASE_DIR}"
echo "ANNOVAR_BASE_DIR: ${ANNOVAR_BASE_DIR}"
echo "OUT_BASE_DIR    : ${OUT_BASE_DIR}"
echo "PLINK2          : ${PLINK2}"
echo "THREADS         : ${THREADS}"
echo

# Determine ancestries to process:
# - If ANCESTRIES env var set, use it.
# - Else auto-detect from ANNOVAR_BASE_DIR (only ancestries with annotation outputs).
declare -a ANC_ARR=()
if [[ -n "${ANCESTRIES:-}" ]]; then
  read -r -a ANC_ARR <<< "${ANCESTRIES}"
else
  mapfile -t ANC_ARR < <(
    find "${ANNOVAR_BASE_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort -u
  )
fi

if [[ "${#ANC_ARR[@]}" -eq 0 ]]; then
  echo "ERROR: No ancestry directories found under: ${ANNOVAR_BASE_DIR}" >&2
  exit 1
fi

echo "ANCESTRIES to process: ${ANC_ARR[*]}"
echo

# Temp files cleanup
TMPFILES=()
cleanup() {
  for f in "${TMPFILES[@]:-}"; do
    [[ -f "$f" ]] && rm -f "$f"
  done
}
trap cleanup EXIT

for ANC in "${ANC_ARR[@]}"; do
  echo "==== ${ANC} ===="

  IN_PFILE_DIR="${PFILE_BASE_DIR}/${ANC}"
  IN_BASE="${IN_PFILE_DIR}/chr11_${ANC}_r10_informative_SORL1"

  if [[ ! -f "${IN_BASE}.pgen" || ! -f "${IN_BASE}.pvar" || ! -f "${IN_BASE}.psam" ]]; then
    echo "  ! Skipping ${ANC}: no SORL1 informative pfile at ${IN_BASE}.{pgen,pvar,psam}"
    echo
    continue
  fi

  ID_LIST="${ANNOVAR_BASE_DIR}/${ANC}/chr11_${ANC}_r10_SORL1_rare_damaging_ids.txt"
  if [[ ! -f "${ID_LIST}" ]]; then
    echo "  ! Skipping ${ANC}: no rare/damaging ID list at ${ID_LIST}"
    echo
    continue
  fi

  # Intersect requested IDs with actual pvar IDs (prevents silent mismatches)
  TMP_REQ="$(mktemp)"
  TMP_PVAR="$(mktemp)"
  TMP_USE="$(mktemp)"
  TMPFILES+=("${TMP_REQ}" "${TMP_PVAR}" "${TMP_USE}")

  # Requested IDs (unique, non-empty)
  awk 'NF{print $1}' "${ID_LIST}" | sort -u > "${TMP_REQ}"

  # Available IDs in pvar (column 3 = ID; skip header)
  awk 'NR>1{print $3}' "${IN_BASE}.pvar" | sort -u > "${TMP_PVAR}"

  # Intersection
  grep -F -x -f "${TMP_REQ}" "${TMP_PVAR}" > "${TMP_USE}" || true

  N_REQ=$(wc -l < "${TMP_REQ}" | awk '{print $1}')
  N_USE=$(wc -l < "${TMP_USE}" | awk '{print $1}')
  N_DROP=$(( N_REQ - N_USE ))

  echo "  * Requested IDs: ${N_REQ}"
  echo "  * Found in pvar : ${N_USE}"
  if [[ "${N_DROP}" -gt 0 ]]; then
    echo "  * Dropped (not in pvar): ${N_DROP}  (expected if a few variants could not be mapped)"
  fi

  if [[ "${N_USE}" -eq 0 ]]; then
    echo "  ! No extractable IDs remain for ${ANC}; skipping."
    echo
    continue
  fi

  # Output ancestry-specific directory
  ANC_OUT_DIR="${OUT_BASE_DIR}/${ANC}"
  mkdir -p "${ANC_OUT_DIR}"
  OUT_BASE="${ANC_OUT_DIR}/chr11_${ANC}_r10_SORL1_rare_damaging"

  echo "  * Extracting rare/damaging SORL1 variants -> ${OUT_BASE}.{pgen,pvar,psam}"
  set +e
  "${PLINK2}" \
    --pfile "${IN_BASE}" \
    --extract "${TMP_USE}" \
    --threads "${THREADS}" \
    --make-pgen \
    --out "${OUT_BASE}"
  rc=$?
  set -e

  if [[ $rc -ne 0 ]]; then
    echo "  ! PLINK2 extraction failed for ${ANC} (exit=${rc}); continuing."
    echo
    continue
  fi

  echo "  * Done ${ANC}."
  echo
done

echo "All requested ancestries processed (rare/damaging SORL1 pfiles)."

