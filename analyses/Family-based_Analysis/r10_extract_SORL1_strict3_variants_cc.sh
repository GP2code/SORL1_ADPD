#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Script-local defaults
# -----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}/cc_data/plink_SORL1_strict3_variants}"
CCDIR="${CCDIR:-${SCRIPT_DIR}/cc_data}"

PLINK2="${PLINK2:-/home/jupyter/bernabe/programs/plink2}"
PLINKBASE="${PLINKBASE:-/home/jupyter/workspace/gp2_tier2_eu_release10/wgs/deepvariant_joint_calling/plink}"

KEEP="${KEEP:-${CCDIR}/r10_wgs_cc_mono_plus_brainbank_UNRELATED_plink_keep.tsv}"

THREADS="${THREADS:-16}"
MEMORY_MB="${MEMORY_MB:-60000}"

mkdir -p "${OUT_DIR}"

echo "OUT_DIR   : ${OUT_DIR}"
echo "CCDIR     : ${CCDIR}"
echo "PLINK2    : ${PLINK2}"
echo "PLINKBASE : ${PLINKBASE}"
echo "KEEP      : ${KEEP}"
echo "THREADS   : ${THREADS}"
echo "MEMORY_MB : ${MEMORY_MB}"
echo

# -----------------------------
# Preconditions
# -----------------------------
if [[ ! -x "${PLINK2}" ]]; then
  echo "ERROR: PLINK2 not found/executable at: ${PLINK2}" >&2
  exit 1
fi

if [[ ! -f "${KEEP}" ]]; then
  echo "ERROR: keep file not found: ${KEEP}" >&2
  exit 1
fi

# Sanity check: keep file should have no header
firstline="$(head -n 1 "${KEEP}" || true)"
if echo "${firstline}" | grep -qiE '(^FID|IID|GP2ID)'; then
  echo "ERROR: keep file appears to have a header. It must be headerless for --keep." >&2
  echo "First line was: ${firstline}" >&2
  exit 1
fi

# -----------------------------
# Variant lists
# -----------------------------
# NOTE: We keep the strict set you decided:
#  - EAS: chr11:121478242:G:A
#  - EUR: chr11:121514222:A:C and chr11:121545392:G:A
# and we are intentionally NOT including the partial-penetrance EAS variant.
EAS_LIST="${OUT_DIR}/EAS_r10_SORL1_strict3_variants.txt"
EUR_LIST="${OUT_DIR}/EUR_r10_SORL1_strict3_variants.txt"

cat > "${EAS_LIST}" << 'EOF'
chr11:121478242:G:A
EOF

cat > "${EUR_LIST}" << 'EOF'
chr11:121514222:A:C
chr11:121545392:G:A
EOF

echo "Wrote variant lists:"
echo "  ${EAS_LIST}:"
cat "${EAS_LIST}" | sed 's/^/    /'
echo "  ${EUR_LIST}:"
cat "${EUR_LIST}" | sed 's/^/    /'
echo

run_extract () {
  local ANC="$1"
  local LIST="$2"
  local INBASE="${PLINKBASE}/${ANC}/chr11_${ANC}_release10"
  local OUTBASE="${OUT_DIR}/chr11_${ANC}_SORL1_strict3_unrelated"

  if [[ ! -f "${INBASE}.pgen" || ! -f "${INBASE}.pvar" || ! -f "${INBASE}.psam" ]]; then
    echo "[${ANC}] WARNING: Missing pfile set at ${INBASE}.{pgen,pvar,psam}. Skipping."
    return 0
  fi

  if [[ ! -s "${LIST}" ]]; then
    echo "[${ANC}] WARNING: Variant list is empty: ${LIST}. Skipping."
    return 0
  fi

  echo "== [${ANC}] Extracting strict variants =="
  "${PLINK2}" \
    --pfile "${INBASE}" \
    --keep "${KEEP}" \
    --extract "${LIST}" \
    --export A \
    --threads "${THREADS}" \
    --memory "${MEMORY_MB}" \
    --out "${OUTBASE}"

  echo "[${ANC}] Wrote:"
  echo "  ${OUTBASE}.raw"
  echo "  ${OUTBASE}.log"
  echo
}

run_extract "EAS" "${EAS_LIST}"
run_extract "EUR" "${EUR_LIST}"

echo "DONE. Outputs in: ${OUT_DIR}"

