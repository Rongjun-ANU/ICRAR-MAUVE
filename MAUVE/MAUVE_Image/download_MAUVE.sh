#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# 1) list of galaxies
targets=(
  NGC4064  NGC4189  NGC4192  NGC4216  NGC4222  NGC4254
  NGC4293  NGC4294  NGC4298  NGC4302  NGC4321  NGC4330
  NGC4351  NGC4380  NGC4383  NGC4388  NGC4394  NGC4396
  NGC4405  NGC4402  NGC4419  NGC4424  NGC4450  IC3392
  NGC4457  NGC4501  NGC4522  NGC4535  NGC4548  NGC4567
  NGC4568  NGC4569  NGC4579  NGC4580  NGC4606  NGC4607
  NGC4654  NGC4689  NGC4694  NGC4698
  M87
)

# 2) resolve name → RA/Dec
get_coords() {
  local name="$1"
  local url
  url=$(printf 'http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/-oI?%s' "$name")
  wget -qO- "$url" \
    | awk '/^%J/ { print $2, $3; exit }'
}

for name in "${targets[@]}"; do
  echo "Resolving $name..."
  coords="$( get_coords "$name" )"
  ra="${coords%% *}"
  dec="${coords##* }"

  echo "Downloading ${name} at RA=${ra}, DEC=${dec}"
  outfile="${name}_Legacy.jpg"
  wget -qO "$outfile" \
    "https://www.legacysurvey.org/viewer/jpeg-cutout?ra=${ra}&dec=${dec}&layer=ls-dr9&pixscale=2.62&width=300&height=300"
  echo "Downloaded $outfile"

  # 3) annotate top-left in red
  # mogrify \
  #   -gravity NorthWest \
  #   -pointsize 24 \
  #   -fill yellow \
  #   -annotate +12+48 "$name" \
  #   "$outfile"
  # echo "Annotated $outfile with label '$name'"
done

echo "All done!"
