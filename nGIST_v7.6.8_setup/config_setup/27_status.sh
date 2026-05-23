#!/usr/bin/env bash
set -uo pipefail

cd "$(dirname "$0")"

RUN_DIR="$(pwd)"
PRODUCT_BASE="/scratch/pawsey1308/mauve/products/v3tk_v7.6.8"
COMPLETION_STRING="MainPipeline: nGIST completed successfully."
STATUS_LOG="27_status_log_$(date +%Y%m%d_%H%M%S).txt"

exec > >(tee "$STATUS_LOG") 2>&1

GALIDS=(
  IC3392
  NGC4064
  NGC4189
  NGC4192
  NGC4293
  NGC4294
  NGC4298
  NGC4302
  NGC4330
  NGC4351
  NGC4383
  NGC4388
  NGC4394
  NGC4396
  NGC4402
  NGC4405
  NGC4419
  NGC4457
  NGC4501
  NGC4522
  NGC4567_8
  NGC4580
  NGC4606
  NGC4607
  NGC4694
  NGC4698
)

timestamp_to_epoch() {
  local stamp="$1"
  if [[ -z "$stamp" ]]; then
    return 1
  fi

  date -d "$stamp" +%s 2>/dev/null
}

format_seconds() {
  local seconds="$1"
  if [[ -z "$seconds" || "$seconds" == "NA" ]]; then
    echo "NA"
    return
  fi

  local days hours minutes
  days=$((seconds / 86400))
  hours=$(((seconds % 86400) / 3600))
  minutes=$(((seconds % 3600) / 60))

  if (( days > 0 )); then
    printf "%dd %02dh %02dm" "$days" "$hours" "$minutes"
  else
    printf "%02dh %02dm" "$hours" "$minutes"
  fi
}

first_timestamp() {
  local log_file="$1"
  awk '/^[0-9][0-9]\/[0-9][0-9]\/[0-9][0-9] [0-9][0-9]:[0-9][0-9]:[0-9][0-9]/ {print substr($0, 1, 17); exit}' "$log_file" 2>/dev/null
}

last_timestamp_for_pattern() {
  local log_file="$1"
  local pattern="$2"
  awk -v pattern="$pattern" '
    $0 ~ pattern && $0 ~ /^[0-9][0-9]\/[0-9][0-9]\/[0-9][0-9] [0-9][0-9]:[0-9][0-9]:[0-9][0-9]/ {
      stamp = substr($0, 1, 17)
    }
    END {
      if (stamp != "") print stamp
    }
  ' "$log_file" 2>/dev/null
}

contains_completion() {
  local log_file="$1"
  [[ -f "$log_file" ]] && tail -n 80 "$log_file" | grep -Fq "$COMPLETION_STRING"
}

extract_spectra_count() {
  local log_file="$1"
  grep -E "Read a total of [0-9]+ spectra" "$log_file" 2>/dev/null | tail -1 | sed -E 's/.*Read a total of ([0-9]+) spectra.*/\1/'
}

extract_bin_count() {
  local log_file="$1"
  grep -E "[0-9]+ Voronoi bins generated" "$log_file" 2>/dev/null | tail -1 | sed -E 's/.* ([0-9]+) Voronoi bins generated.*/\1/'
}

extract_gas_work_units() {
  local log_file="$1"
  local bin_count spectra_count

  bin_count="$(extract_bin_count "$log_file")"
  spectra_count="$(extract_spectra_count "$log_file")"

  if [[ "$bin_count" =~ ^[0-9]+$ && "$spectra_count" =~ ^[0-9]+$ ]]; then
    echo $((bin_count + spectra_count))
  elif [[ "$spectra_count" =~ ^[0-9]+$ ]]; then
    echo "$spectra_count"
  elif [[ "$bin_count" =~ ^[0-9]+$ ]]; then
    echo "$bin_count"
  fi
}

last_info_line() {
  local log_file="$1"
  grep -E " - INFO[[:space:]]+- " "$log_file" 2>/dev/null | tail -1 | sed -E 's/[[:space:]]+/ /g'
}

job_state() {
  local galid="$1"
  local job_name="${galid}_v3tk_v7.6.8"
  local queue_user="${USER:-${LOGNAME:-}}"

  if ! command -v squeue >/dev/null 2>&1; then
    return 0
  fi

  if [[ -n "$queue_user" ]]; then
    squeue -h -u "$queue_user" -o "%j|%T|%i|%M|%L" 2>/dev/null
  else
    squeue -h -o "%j|%T|%i|%M|%L" 2>/dev/null
  fi | awk -F'|' -v job_name="$job_name" '
    $1 == job_name {
      print $2 " jobid=" $3 " elapsed=" $4 " left=" $5
      exit
    }
  '
}

classify_status() {
  local galid="$1"
  local product_log="${PRODUCT_BASE}/${galid}/LOGFILE"
  local run_log="${RUN_DIR}/${galid}_v3tk_v7.6.8.log"
  local state

  if contains_completion "$product_log"; then
    echo "FINISHED"
    return
  fi

  state="$(job_state "$galid")"
  if [[ -n "$state" ]]; then
    echo "RUNNING"
    return
  fi

  if [[ -f "$run_log" ]] && grep -Fq "DUE TO TIME LIMIT" "$run_log"; then
    echo "TIMEOUT_RESBATCH"
    return
  fi

  if [[ -f "$run_log" && ! -s "$run_log" ]]; then
    echo "RUNNING_EMPTY_LOG"
    return
  fi

  if [[ -f "$product_log" ]]; then
    echo "UNFINISHED_UNKNOWN"
    return
  fi

  echo "MISSING_LOGS"
}

stage_label() {
  local log_file="$1"

  if [[ ! -f "$log_file" ]]; then
    echo "no product LOGFILE"
  elif contains_completion "$log_file"; then
    echo "completed"
  elif grep -Eiq "_emission|emissionLines|Produced .*gas|gas.*maps" "$log_file"; then
    echo "gas/emission stage reached"
  elif grep -Eq "_continuumCube|ppxf_cont_wrapper" "$log_file"; then
    echo "continuum stage reached"
  elif grep -Fq "Produced stellar kinematics maps" "$log_file"; then
    echo "stellar kinematics maps done"
  elif grep -Fq "_stellarKinematics" "$log_file"; then
    echo "stellar kinematics stage reached"
  elif grep -Fq "_bin_spectra.hdf5" "$log_file"; then
    echo "binned log spectra written"
  elif grep -Fq "_all_spectra.hdf5" "$log_file"; then
    echo "all log spectra written"
  elif grep -Fq "Voronoi bins generated" "$log_file"; then
    echo "Voronoi binning done"
  elif grep -Fq "Finished reading the MUSE cube" "$log_file"; then
    echo "cube read done"
  else
    echo "started or unknown"
  fi
}

stage_pattern() {
  local log_file="$1"

  if [[ ! -f "$log_file" ]]; then
    echo ""
  elif grep -Eiq "_emission|emissionLines|Produced .*gas|gas.*maps" "$log_file"; then
    echo "_emission|emissionLines|Produced .*gas|gas.*maps"
  elif grep -Eq "_continuumCube|ppxf_cont_wrapper" "$log_file"; then
    echo "_continuumCube|ppxf_cont_wrapper"
  elif grep -Fq "Produced stellar kinematics maps" "$log_file"; then
    echo "Produced stellar kinematics maps"
  elif grep -Fq "_stellarKinematics" "$log_file"; then
    echo "_stellarKinematics"
  elif grep -Fq "_bin_spectra.hdf5" "$log_file"; then
    echo "_bin_spectra.hdf5"
  elif grep -Fq "_all_spectra.hdf5" "$log_file"; then
    echo "_all_spectra.hdf5"
  elif grep -Fq "Voronoi bins generated" "$log_file"; then
    echo "Voronoi bins generated"
  elif grep -Fq "Finished reading the MUSE cube" "$log_file"; then
    echo "Finished reading the MUSE cube"
  else
    echo ""
  fi
}

max_from_stdin() {
  sort -n | awk '
    NF {
      max = $1
      seen = 1
    }
    END {
      if (!seen) {
        exit
      }
      print int(max + 0.5)
    }
  '
}

estimate_remaining_seconds() {
  local pattern="$1"
  local target_work="$2"
  local completion_ts marker_ts completion_epoch marker_epoch duration ref_work finished_log

  if [[ -z "$pattern" || ${#FINISHED_LOGS[@]} -eq 0 ]]; then
    echo ""
    return
  fi

  for finished_log in "${FINISHED_LOGS[@]}"; do
    marker_ts="$(last_timestamp_for_pattern "$finished_log" "$pattern")"
    completion_ts="$(last_timestamp_for_pattern "$finished_log" "$COMPLETION_STRING")"

    if [[ -z "$marker_ts" || -z "$completion_ts" ]]; then
      continue
    fi

    marker_epoch="$(timestamp_to_epoch "$marker_ts" || true)"
    completion_epoch="$(timestamp_to_epoch "$completion_ts" || true)"

    if [[ -z "$marker_epoch" || -z "$completion_epoch" ]]; then
      continue
    fi

    duration=$((completion_epoch - marker_epoch))
    if (( duration <= 0 )); then
      continue
    fi

    ref_work="$(extract_gas_work_units "$finished_log")"
    if [[ "$target_work" =~ ^[0-9]+$ && "$ref_work" =~ ^[0-9]+$ && "$ref_work" -gt 0 ]]; then
      awk -v duration="$duration" -v target_work="$target_work" -v ref_work="$ref_work" 'BEGIN {print int(duration * target_work / ref_work + 0.5)}'
    else
      echo "$duration"
    fi
  done | max_from_stdin
}

FINISHED_LOGS=()
for galid in "${GALIDS[@]}"; do
  product_log="${PRODUCT_BASE}/${galid}/LOGFILE"
  if contains_completion "$product_log"; then
    FINISHED_LOGS+=("$product_log")
  fi
done

finished_count=0
running_count=0
timeout_count=0
unknown_count=0
missing_count=0
resbatch_list=()

echo "nGIST v7.6.8 batch status"
echo "Generated: $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "Run directory: ${RUN_DIR}"
echo "Product directory: ${PRODUCT_BASE}"
echo "Status log: ${RUN_DIR}/${STATUS_LOG}"
echo

printf "%-12s %-17s %-34s %-17s %-9s %-8s %-10s %-14s %s\n" \
  "GALID" "STATUS" "SQUEUE" "LAST_TIME" "SPECTRA" "BINS" "GAS_WORK" "EST_REMAIN" "STAGE_OR_ACTION"
printf "%-12s %-17s %-34s %-17s %-9s %-8s %-10s %-14s %s\n" \
  "------------" "-----------------" "----------------------------------" "-----------------" "---------" "--------" "----------" "--------------" "----------------"

for galid in "${GALIDS[@]}"; do
  product_log="${PRODUCT_BASE}/${galid}/LOGFILE"
  run_log="${RUN_DIR}/${galid}_v3tk_v7.6.8.log"
  status="$(classify_status "$galid")"
  state="$(job_state "$galid")"
  stage="$(stage_label "$product_log")"
  spectra="$(extract_spectra_count "$product_log")"
  bins="$(extract_bin_count "$product_log")"
  gas_work="$(extract_gas_work_units "$product_log")"
  last_ts="$(last_timestamp_for_pattern "$product_log" ".*")"
  estimate="NA"
  action="$stage"

  [[ -z "$state" ]] && state="-"
  [[ -z "$spectra" ]] && spectra="-"
  [[ -z "$bins" ]] && bins="-"
  [[ -z "$gas_work" ]] && gas_work="-"
  [[ -z "$last_ts" ]] && last_ts="-"

  case "$status" in
    FINISHED)
      finished_count=$((finished_count + 1))
      action="done"
      ;;
    RUNNING|RUNNING_EMPTY_LOG)
      running_count=$((running_count + 1))
      if [[ "$status" == "RUNNING_EMPTY_LOG" ]]; then
        action="not complete; run log is empty, likely still running or just started"
      else
        action="not complete; still in scheduler"
      fi
      ;;
    TIMEOUT_RESBATCH)
      timeout_count=$((timeout_count + 1))
      resbatch_list+=("$galid")
      pattern="$(stage_pattern "$product_log")"
      estimate_seconds="$(estimate_remaining_seconds "$pattern" "$gas_work")"
      if [[ -n "$estimate_seconds" ]]; then
        estimate="$(format_seconds "$estimate_seconds")"
      fi
      action="timeout; resubmit with sbatch ${galid}_v3tk_v7.6.8_setonix.slurm; latest stage: ${stage}"
      ;;
    MISSING_LOGS)
      missing_count=$((missing_count + 1))
      action="missing product LOGFILE and run log context"
      ;;
    *)
      unknown_count=$((unknown_count + 1))
      action="not complete; inspect product LOGFILE and ${run_log}"
      ;;
  esac

  printf "%-12s %-17s %-34s %-17s %-9s %-8s %-10s %-14s %s\n" \
    "$galid" "$status" "$state" "$last_ts" "$spectra" "$bins" "$gas_work" "$estimate" "$action"
done

echo
echo "Summary"
echo "  Finished: ${finished_count}"
echo "  Still running or queued: ${running_count}"
echo "  Timed out and should be resubmitted: ${timeout_count}"
echo "  Unfinished with unknown state: ${unknown_count}"
echo "  Missing logs: ${missing_count}"

if (( timeout_count > 0 )); then
  echo
  echo "Resubmit commands for timeout jobs"
  for galid in "${resbatch_list[@]}"; do
    echo "  sbatch ${galid}_v3tk_v7.6.8_setonix.slurm"
  done
fi

echo
echo "Notes"
echo "  Finished means the product LOGFILE contains: ${COMPLETION_STRING}"
echo "  TIMEOUT_RESBATCH means the run log contains: DUE TO TIME LIMIT"
echo "  GAS_WORK = SPECTRA + BINS, matching GAS LEVEL=BOTH spaxel-level plus bin-level fitting."
echo "  Estimates use the maximum finished-job remainder at the same reached stage, scaled by GAS_WORK when available."
echo "  Because nGIST can skip completed modules on restart, timeout estimates are for the remaining pipeline from the latest detected stage."
