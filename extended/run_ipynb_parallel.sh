#!/usr/bin/env bash
set -euo pipefail

# --- 0. Start timing ----------------------------------------------------------
start_time=$(date +%s)

# --- 1. Activate Conda -------------------------------------------------------
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ICRAR          # your environment name

# --- 2. Runtime parameters ---------------------------------------------------
export THREADS="${THREADS:-$(sysctl -n hw.ncpu)}"   # default = all CPU cores
export TIMEOUT=-1                                   # disable cell timeout

# --- 3. Parallel execution ----------------------------------------------------
parallel --will-cite --bar -j "$THREADS" --halt now,fail=1 '
  nb="{}"
  log="${nb%.ipynb}.log"

  echo "▶ Executing $nb" >&2

  # ---- run notebook and capture exit status ---------------------------------
  jupyter nbconvert --execute --inplace \
      --ExecutePreprocessor.timeout=$TIMEOUT "$nb" \
      >"$log" 2>&1
  status=$?

  # ---- if successful, delete its log; else keep it for debugging ------------
  if [[ $status -eq 0 ]]; then
      rm -f "$log"
  fi

  exit $status               # propagate status back to GNU parallel
' ::: Line_to_SFR_*.ipynb

# --- 4. Final summary ---------------------------------------------------------
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "✅ All notebooks executed successfully."
echo "⏱️  Total runtime: ${runtime} seconds"
