#!/usr/bin/env bash
# 2026-10-30: Short entry point for the existing SFH16 Mass -> QC loop.
# Full sample: 40 galaxies in 39 cubes (NGC4567_8 is one combined cube).
# Run from further with the science Python environment active on CANFAR.
#
# Usage:
#   bash 7000_SFH16.sh --dry-run       # preview the full selection
#   bash 7000_SFH16.sh                 # process all available SFH16 runs
#   bash 7000_SFH16.sh NGC4689         # process selected runs only
#   bash 7000_SFH16.sh --help          # all options
# Implementation and detailed instructions: QC_ngist_v3tk_v768_7000_SFH16.{sh,md}
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$script_dir/QC_ngist_v3tk_v768_7000_SFH16.sh" "$@"
