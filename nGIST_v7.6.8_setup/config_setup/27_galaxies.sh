#!/usr/bin/env bash

ALL_GALIDS=(
  IC3392
  NGC4064
  NGC4189
  NGC4192
  NGC4216
  NGC4222
  NGC4254
  NGC4293
  NGC4294
  NGC4298
  NGC4302
  NGC4321
  NGC4330
  NGC4351
  NGC4380
  NGC4383
  NGC4388
  NGC4394
  NGC4396
  NGC4405
  NGC4402
  NGC4419
  NGC4424
  NGC4450
  NGC4457
  NGC4501
  NGC4522
  NGC4535
  NGC4548
  NGC4567_8
  NGC4569
  NGC4579
  NGC4580
  NGC4606
  NGC4607
  NGC4654
  NGC4689
  NGC4694
  NGC4698
)

is_known_galid() {
  local candidate="$1"
  local galid

  for galid in "${ALL_GALIDS[@]}"; do
    if [[ "$candidate" == "$galid" ]]; then
      return 0
    fi
  done

  return 1
}

looks_like_galid() {
  local candidate="$1"
  [[ "$candidate" =~ ^(IC|NGC)[0-9][0-9_]*$ ]]
}

select_galids() {
  SELECTED_GALIDS=()

  if [[ $# -eq 0 ]]; then
    SELECTED_GALIDS=("${ALL_GALIDS[@]}")
    return 0
  fi

  local requested
  for requested in "$@"; do
    if ! is_known_galid "$requested"; then
      echo "ERROR: unknown galaxy ID: $requested" >&2
      echo "Known galaxy IDs: ${ALL_GALIDS[*]}" >&2
      return 2
    fi

    SELECTED_GALIDS+=("$requested")
  done
}
