#!/usr/bin/env bash
set -uo pipefail

cd "$(dirname "$0")"

RUN_DIR="$(pwd)"
PRODUCT_BASE_REGULAR="/scratch/pawsey1308/mauve/products/v3tk_v7.6.8"
PRODUCT_BASE_7000="/scratch/pawsey1308/mauve/products/v3tk_v7.6.8_7000"
COMPLETION_STRING="MainPipeline: nGIST completed successfully."
STATUS_LOG="27_status_log_$(date +%Y%m%d_%H%M%S).txt"
LONG_QUEUE_ESTIMATE_WARNING_SECONDS=$((22 * 3600))
RUN_VARIANTS=(
  "regular||${PRODUCT_BASE_REGULAR}"
  "7000|_7000|${PRODUCT_BASE_7000}"
)

. ./27_galaxies.sh
select_galids "$@"
GALIDS=("${SELECTED_GALIDS[@]}")

timestamp_to_epoch() {
  local stamp="$1"
  if [[ -z "$stamp" ]]; then
    return 1
  fi

  date -d "$stamp" +%s 2>/dev/null || date -j -f "%m/%d/%y %H:%M:%S" "$stamp" +%s 2>/dev/null
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

config_method() {
  local config_file="$1"
  local section="$2"

  if [[ ! -f "$config_file" ]]; then
    return
  fi

  awk -v section="$section" '
    $0 ~ "^" section ":" {
      in_section = 1
      next
    }
    in_section && $0 ~ /^[A-Z0-9_]+:/ {
      exit
    }
    in_section && $0 ~ /^[[:space:]]+METHOD:/ {
      value = $0
      sub(/^[[:space:]]+METHOD:[[:space:]]*/, "", value)
      gsub(/[[:space:]]+$/, "", value)
      print value
      exit
    }
  ' "$config_file" 2>/dev/null
}

method_is_enabled() {
  local method="$1"
  local normalised

  normalised="$(printf "%s" "$method" | tr '[:upper:]' '[:lower:]')"
  [[ -n "$normalised" && "$normalised" != "false" && "$normalised" != "none" && "$normalised" != "null" && "$normalised" != "off" ]]
}

module_is_enabled() {
  local config_file="$1"
  local module="$2"
  local method

  method="$(config_method "$config_file" "$module")"
  method_is_enabled "$method"
}

last_info_line() {
  local log_file="$1"
  grep -E " - INFO[[:space:]]+- " "$log_file" 2>/dev/null | tail -1 | sed -E 's/[[:space:]]+/ /g'
}

job_state() {
  local galid="$1"
  local job_name="${galid}_v3tk_v7.6.8${RUN_SUFFIX}"
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

slurm_partition() {
  local galid="$1"
  local slurm_file="${RUN_DIR}/${galid}_v3tk_v7.6.8${RUN_SUFFIX}_setonix.slurm"

  awk -F= '/^#SBATCH --partition=/ {print $2; exit}' "$slurm_file" 2>/dev/null
}

slurm_walltime() {
  local galid="$1"
  local slurm_file="${RUN_DIR}/${galid}_v3tk_v7.6.8${RUN_SUFFIX}_setonix.slurm"

  awk -F= '/^#SBATCH --time=/ {print $2; exit}' "$slurm_file" 2>/dev/null
}

queue_warning_message() {
  local galid="$1"
  local reason="$2"
  local partition walltime

  partition="$(slurm_partition "$galid")"
  walltime="$(slurm_walltime "$galid")"

  if [[ "$partition" == "work" ]]; then
    echo "WARNING: ${reason}; update slurm from partition=work to partition=long time=96:00:00"
  elif [[ "$partition" == "long" || "$partition" == "highmem" ]]; then
    echo "WARNING: ${reason}; current slurm uses partition=${partition} time=${walltime:-unknown}, do not resubmit on work"
  else
    echo "WARNING: ${reason}; use a 96h queue, normally partition=long unless highmem is required"
  fi
}

append_warning_text() {
  local current="$1"
  local extra="$2"

  if [[ -n "$current" ]]; then
    printf "%s; %s" "$current" "$extra"
  else
    printf "%s" "$extra"
  fi
}

classify_status() {
  local galid="$1"
  local product_log="${PRODUCT_BASE}/${galid}/LOGFILE"
  local run_log="${RUN_DIR}/${galid}_v3tk_v7.6.8${RUN_SUFFIX}.log"
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
  elif grep -Fq "Produced SFH maps" "$log_file"; then
    echo "SFH maps done"
  elif grep -Eq "_starFormationHistories|ppxf_sfh_wrapper" "$log_file"; then
    echo "SFH stage reached"
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

count_pattern() {
  local log_file="$1"
  local pattern="$2"

  awk -v pattern="$pattern" '$0 ~ pattern {count++} END {print count + 0}' "$log_file" 2>/dev/null
}

last_line_number_for_pattern() {
  local log_file="$1"
  local pattern="$2"

  awk -v pattern="$pattern" '$0 ~ pattern {line = NR} END {if (line != "") print line}' "$log_file" 2>/dev/null
}

has_pattern_after_line() {
  local log_file="$1"
  local line_number="$2"
  local pattern="$3"

  awk -v line_number="$line_number" -v pattern="$pattern" '
    NR > line_number && $0 ~ pattern {
      found = 1
      exit
    }
    END {
      exit(found ? 0 : 1)
    }
  ' "$log_file" 2>/dev/null
}

repeated_unfinished_module_label() {
  local log_file="$1"
  local count last_start

  if [[ ! -f "$log_file" ]] || contains_completion "$log_file"; then
    return
  fi

  count="$(count_pattern "$log_file" "_emissionLines: Using the emissionLines routine")"
  if (( count >= 2 )); then
    last_start="$(last_line_number_for_pattern "$log_file" "_emissionLines: Using the emissionLines routine")"
    if [[ -n "$last_start" ]] && ! has_pattern_after_line "$log_file" "$last_start" "Emission Line Fitting done|_starFormationHistories: Using the starFormationHistories routine|MainPipeline: nGIST completed successfully"; then
      echo "emissionLines"
      return
    fi
  fi

  count="$(count_pattern "$log_file" "_starFormationHistories: Using the starFormationHistories routine")"
  if (( count >= 2 )); then
    last_start="$(last_line_number_for_pattern "$log_file" "_starFormationHistories: Using the starFormationHistories routine")"
    if [[ -n "$last_start" ]] && ! has_pattern_after_line "$log_file" "$last_start" "Produced SFH maps|_lineStrengths:|_userModules:|MainPipeline: nGIST completed successfully"; then
      echo "starFormationHistories"
      return
    fi
  fi
}

module_is_complete() {
  local module="$1"
  local log_file="$2"

  case "$module" in
    GAS)
      grep -Fq "Emission Line Fitting done" "$log_file" 2>/dev/null || contains_completion "$log_file"
      ;;
    SFH)
      grep -Fq "Produced SFH maps" "$log_file" 2>/dev/null || contains_completion "$log_file"
      ;;
    LS)
      grep -Fq "Produced line-strength maps" "$log_file" 2>/dev/null || contains_completion "$log_file"
      ;;
    UMOD)
      contains_completion "$log_file"
      ;;
    *)
      return 1
      ;;
  esac
}

gas_resume_level() {
  local galid="$1"
  local log_file="$2"
  local product_dir="${PRODUCT_BASE}/${galid}"

  if [[ -f "${product_dir}/${galid}_gas_bin.fits" ]] || grep -Fq "${galid}_gas_bin.fits" "$log_file" 2>/dev/null; then
    echo "SPAXEL"
  else
    echo "BOTH"
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

estimate_module_seconds() {
  local module="$1"
  local galid="$2"
  local log_file="$3"
  local start_pattern end_pattern target_work ref_work finished_log
  local start_ts end_ts start_epoch end_epoch duration gas_level

  case "$module" in
    GAS)
      gas_level="$(gas_resume_level "$galid" "$log_file")"
      if [[ "$gas_level" == "SPAXEL" ]]; then
        start_pattern="Using full spectral library for ppxf on SPAXEL level"
        target_work="$(extract_spectra_count "$log_file")"
      else
        start_pattern="_emissionLines: Using the emissionLines routine"
        target_work="$(extract_gas_work_units "$log_file")"
      fi
      end_pattern="Emission Line Fitting done"
      ;;
    SFH)
      start_pattern="_starFormationHistories: Using the starFormationHistories routine"
      end_pattern="Produced SFH maps"
      target_work="$(extract_bin_count "$log_file")"
      ;;
    *)
      echo ""
      return
      ;;
  esac

  if [[ -z "$start_pattern" || -z "$end_pattern" || ${#FINISHED_LOGS[@]} -eq 0 ]]; then
    echo ""
    return
  fi

  for finished_log in "${FINISHED_LOGS[@]}"; do
    start_ts="$(last_timestamp_for_pattern "$finished_log" "$start_pattern")"
    end_ts="$(last_timestamp_for_pattern "$finished_log" "$end_pattern")"

    if [[ -z "$start_ts" || -z "$end_ts" ]]; then
      continue
    fi

    start_epoch="$(timestamp_to_epoch "$start_ts" || true)"
    end_epoch="$(timestamp_to_epoch "$end_ts" || true)"

    if [[ -z "$start_epoch" || -z "$end_epoch" ]]; then
      continue
    fi

    duration=$((end_epoch - start_epoch))
    if (( duration <= 0 )); then
      continue
    fi

    case "$module" in
      GAS)
        if [[ "$(gas_resume_level "$galid" "$log_file")" == "SPAXEL" ]]; then
          ref_work="$(extract_spectra_count "$finished_log")"
        else
          ref_work="$(extract_gas_work_units "$finished_log")"
        fi
        ;;
      SFH)
        ref_work="$(extract_bin_count "$finished_log")"
        ;;
    esac

    if [[ "$target_work" =~ ^[0-9]+$ && "$ref_work" =~ ^[0-9]+$ && "$ref_work" -gt 0 ]]; then
      awk -v duration="$duration" -v target_work="$target_work" -v ref_work="$ref_work" 'BEGIN {print int(duration * target_work / ref_work + 0.5)}'
    else
      echo "$duration"
    fi
  done | max_from_stdin
}

remaining_module_estimates() {
  local galid="$1"
  local log_file="$2"
  local config_file="$3"
  local module module_seconds module_text
  local total_seconds=0
  local has_estimate=0
  local has_remaining=0
  local over_limit=0
  local breakdown=""

  if [[ ! -f "$log_file" || ! -f "$config_file" ]] || contains_completion "$log_file"; then
    return
  fi

  for module in GAS SFH LS UMOD; do
    if ! module_is_enabled "$config_file" "$module"; then
      continue
    fi

    if module_is_complete "$module" "$log_file"; then
      continue
    fi

    has_remaining=1
    module_seconds="$(estimate_module_seconds "$module" "$galid" "$log_file")"

    if [[ "$module_seconds" =~ ^[0-9]+$ ]]; then
      has_estimate=1
      total_seconds=$((total_seconds + module_seconds))
      module_text="${module}~$(format_seconds "$module_seconds")"
      if (( module_seconds > LONG_QUEUE_ESTIMATE_WARNING_SECONDS )); then
        over_limit=1
      fi
    else
      module_text="${module}~NA"
    fi

    if [[ -n "$breakdown" ]]; then
      breakdown="${breakdown} + ${module_text}"
    else
      breakdown="$module_text"
    fi
  done

  if (( has_remaining == 0 )); then
    return
  fi

  if (( has_estimate == 1 )); then
    printf "%s|%s|%s\n" "$total_seconds" "$breakdown" "$over_limit"
  else
    printf "|%s|%s\n" "$breakdown" "$over_limit"
  fi
}

{
FINISHED_LOGS=()
for run_variant in "${RUN_VARIANTS[@]}"; do
  IFS='|' read -r RUN_LABEL RUN_SUFFIX PRODUCT_BASE <<< "$run_variant"
  for galid in "${GALIDS[@]}"; do
    product_log="${PRODUCT_BASE}/${galid}/LOGFILE"
    if contains_completion "$product_log"; then
      FINISHED_LOGS+=("$product_log")
    fi
  done
done

finished_count=0
running_count=0
timeout_count=0
unknown_count=0
missing_count=0
resbatch_list=()
long_queue_warning_list=()

echo "nGIST v7.6.8 batch status"
echo "Generated: $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "Run directory: ${RUN_DIR}"
echo "Product directories:"
echo "  regular: ${PRODUCT_BASE_REGULAR}"
echo "  7000: ${PRODUCT_BASE_7000}"
echo "Status log: ${RUN_DIR}/${STATUS_LOG}"
echo

printf "%-8s %-12s %-17s %-34s %-17s %-9s %-8s %-10s %-14s %s\n" \
  "RUN" "GALID" "STATUS" "SQUEUE" "LAST_TIME" "SPECTRA" "BINS" "GAS_WORK" "EST_REMAIN" "STAGE_OR_ACTION"
printf "%-8s %-12s %-17s %-34s %-17s %-9s %-8s %-10s %-14s %s\n" \
  "--------" "------------" "-----------------" "----------------------------------" "-----------------" "---------" "--------" "----------" "--------------" "----------------"

for run_variant in "${RUN_VARIANTS[@]}"; do
  IFS='|' read -r RUN_LABEL RUN_SUFFIX PRODUCT_BASE <<< "$run_variant"
  for galid in "${GALIDS[@]}"; do
    product_log="${PRODUCT_BASE}/${galid}/LOGFILE"
    config_file="${PRODUCT_BASE}/${galid}/CONFIG"
    run_log="${RUN_DIR}/${galid}_v3tk_v7.6.8${RUN_SUFFIX}.log"
    status="$(classify_status "$galid")"
    state="$(job_state "$galid")"
    stage="$(stage_label "$product_log")"
    spectra="$(extract_spectra_count "$product_log")"
    bins="$(extract_bin_count "$product_log")"
    gas_work="$(extract_gas_work_units "$product_log")"
    last_ts="$(last_timestamp_for_pattern "$product_log" ".*")"
    estimate="NA"
    estimate_seconds=""
    action="$stage"
    warning_text=""
    repeated_module="$(repeated_unfinished_module_label "$product_log")"
    remaining_info="$(remaining_module_estimates "$galid" "$product_log" "$config_file")"
    remaining_total_seconds=""
    remaining_breakdown=""
    remaining_over_limit=0

    if [[ -n "$remaining_info" ]]; then
      IFS='|' read -r remaining_total_seconds remaining_breakdown remaining_over_limit <<< "$remaining_info"
      if [[ "$remaining_total_seconds" =~ ^[0-9]+$ ]]; then
        estimate_seconds="$remaining_total_seconds"
        estimate="$(format_seconds "$estimate_seconds")"
      fi
    fi

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
        resbatch_list+=("${galid}|${RUN_SUFFIX}")
        action="timeout; resubmit with sbatch ${galid}_v3tk_v7.6.8${RUN_SUFFIX}_setonix.slurm; latest stage: ${stage}"
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

    if [[ -n "$remaining_breakdown" && "$status" != "FINISHED" ]]; then
      action="${action}; remaining active modules: ${remaining_breakdown}"
    fi

    if [[ -n "$repeated_module" ]]; then
      warning_text="$(queue_warning_message "$galid" "repeated unfinished ${repeated_module} module")"
    fi

    if [[ "$remaining_over_limit" == "1" ]]; then
      warning_text="$(append_warning_text "$warning_text" "$(queue_warning_message "$galid" "at least one remaining module estimate > 22h")")"
    elif [[ "$estimate_seconds" =~ ^[0-9]+$ ]] && (( estimate_seconds > LONG_QUEUE_ESTIMATE_WARNING_SECONDS )); then
      warning_text="$(append_warning_text "$warning_text" "$(queue_warning_message "$galid" "EST_REMAIN > 22h")")"
    fi

    if [[ -n "$warning_text" ]]; then
      action="${action}; ${warning_text}"
      long_queue_warning_list+=("${RUN_LABEL} ${galid}: ${warning_text}")
    fi

    printf "%-8s %-12s %-17s %-34s %-17s %-9s %-8s %-10s %-14s %s\n" \
      "$RUN_LABEL" "$galid" "$status" "$state" "$last_ts" "$spectra" "$bins" "$gas_work" "$estimate" "$action"
  done
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
  for resbatch_item in "${resbatch_list[@]}"; do
    IFS='|' read -r galid run_suffix <<< "$resbatch_item"
    echo "  sbatch ${galid}_v3tk_v7.6.8${run_suffix}_setonix.slurm"
  done
fi

if (( ${#long_queue_warning_list[@]} > 0 )); then
  echo
  echo "Long-queue warnings"
  for warning in "${long_queue_warning_list[@]}"; do
    echo "  ${warning}"
  done
fi

echo
echo "Notes"
echo "  Finished means the product LOGFILE contains: ${COMPLETION_STRING}"
echo "  TIMEOUT_RESBATCH means the run log contains: DUE TO TIME LIMIT"
echo "  GAS_WORK = SPECTRA + BINS, matching GAS LEVEL=BOTH spaxel-level plus bin-level fitting."
echo "  EST_REMAIN sums enabled, unfinished modules from each product CONFIG."
echo "  GAS estimates scale by SPECTRA for SPAXEL-only resume, or SPECTRA+BINS before BIN gas is complete."
echo "  SFH estimates scale by BINS. LS and UMOD are listed as NA if enabled but no estimator is implemented."
echo "  Estimates use the maximum scaled module time from comparable finished jobs."
echo "  Long-queue warnings mean the same module restarted without completing, one module estimate is longer than 22h, or total EST_REMAIN is longer than 22h."
} 2>&1 | tee "$STATUS_LOG"
