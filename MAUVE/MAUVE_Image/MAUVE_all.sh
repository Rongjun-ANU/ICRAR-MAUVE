#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# MAUVE_all.sh — 300x300 stamps @ 2.62"/px, PNG + LARGE LABELS
# - Accurate TAN projection @ 2.62"/px with EAST=LEFT (RA flipped)
# - Overlaps: earlier target is masked under later one
# - Labels: placed ONLY on empty background (avoid stamps & labels)
# - Output: lossless PNG with very large, crisp text
# ============================================================

# ---------- Config ----------
SCALE_AS_PER_PX=2.62
STAMP_PX=300
OUTFILE="MAUVE_all.png"     # PNG output (lossless)

# ***** Label geometry (5× larger) *****
LABEL_PT=140                # point size (5× the previous 28)
LABEL_CHAR_W=80             # ~0.57*LABEL_PT; good heuristic for width per char
LABEL_PAD_X=50              # 5× the previous 10
LABEL_PAD_Y=40              # 5× the previous 8
LABEL_H=$(( LABEL_PT + 2*LABEL_PAD_Y ))
# Try bigger search margins first, since labels are large:
LABEL_MARGIN_LIST="40,80,120,160,240,320,480,640,800,1000,1200,1400"

# modest IM limits
export MAGICK_TMPDIR="${PWD}/_magick_tmp"
mkdir -p "$MAGICK_TMPDIR"
IM_LIMITS=( -limit memory 1GiB -limit map 8GiB -limit disk 50GiB )

# ---------- Tools ----------
IM=""
if command -v magick >/dev/null 2>&1; then IM="magick"
elif command -v convert >/dev/null 2>&1; then IM="convert"
else echo "Error: ImageMagick not found."; exit 1; fi
command -v wget >/dev/null 2>&1 || { echo "Error: wget not found."; exit 1; }
command -v awk  >/dev/null 2>&1 || { echo "Error: awk not found.";  exit 1; }

# ---------- Targets (order controls who is on top in overlaps) ----------
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

# ---------- Resolve RA/Dec ----------
get_coords() {
  local name="$1" url
  url=$(printf 'http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/-oI?%s' "$name")
  wget -qO- "$url" | awk '/^%J/ { print $2, $3; exit }'
}

COORDS_FILE="coords.tsv"; : > "$COORDS_FILE"
for name in "${targets[@]}"; do
  coords="$( get_coords "$name" || true )"
  [[ -z "${coords// }" ]] && { echo "Failed to resolve $name" >&2; exit 1; }
  printf "%s\t%s\t%s\n" "$name" "${coords%% *}" "${coords##* }" >> "$COORDS_FILE"
done

# ---------- TAN projection @ fixed 2.62"/px (EAST LEFT: flip RA) ----------
STAMP_ARCSEC=$(awk -v s=$SCALE_AS_PER_PX -v n=$STAMP_PX 'BEGIN{printf "%.6f", s*n}')
MANIFEST="manifest.txt"
awk -v SCALE="$SCALE_AS_PER_PX" -v STAMP_ARCSEC="$STAMP_ARCSEC" -v STAMP_PX="$STAMP_PX" '
BEGIN{ pi=4*atan2(1,1); d2r=pi/180.0; r2as=206264.8062470964 }
{
  name[NR]=$1; ra_deg[NR]=$2; dec_deg[NR]=$3
  ra_rad[NR]=ra_deg[NR]*d2r; dec_rad[NR]=dec_deg[NR]*d2r
  sra+=sin(ra_rad[NR]); cra+=cos(ra_rad[NR]); sdec+=dec_rad[NR]
}
END{
  n=NR; if(n==0) exit 1
  ra0=atan2(sra/n,cra/n); dec0=sdec/n; sd0=sin(dec0); cd0=cos(dec0)

  minx=1e99; maxx=-1e99; miny=1e99; maxy=-1e99
  for(i=1;i<=n;i++){
    dra=ra_rad[i]-ra0; if(dra>pi) dra-=2*pi; if(dra<-pi) dra+=2*pi
    sd=sin(dec_rad[i]); cd=cos(dec_rad[i])
    cosc=sd0*sd + cd0*cd*cos(dra); if(cosc==0) next
    x_as=(cd*sin(dra))/cosc*r2as
    y_as=(cd0*sd - sd0*cd*cos(dra))/cosc*r2as
    if(x_as<minx)minx=x_as; if(x_as>maxx)maxx=x_as
    if(y_as<miny)miny=y_as; if(y_as>maxy)maxy=y_as
    x[i]=x_as; y[i]=y_as
  }
  xspan=maxx-minx; yspan=maxy-miny
  margin=STAMP_ARCSEC/2.0 + 100.0
  totalX=xspan+2*margin; totalY=yspan+2*margin
  W=int(totalX/SCALE + 0.5); H=int(totalY/SCALE + 0.5); if(W<=0||H<=0) exit 1

  printf "CANVAS\t%d\t%d\t%.6f\t%d\n", W, H, SCALE, STAMP_PX > "'$MANIFEST'"

  for(i=1;i<=n;i++){
    # Flip RA horizontally so EAST is to the LEFT (match Legacy stamps)
    px=int(((maxx - x[i] + margin)/SCALE) + 0.5)
    py=int(((maxy - y[i] + margin)/SCALE) + 0.5)  # north up
    printf "%s\t%d\t%d\n", name[i], px, py >> "'$MANIFEST'"
  }
}' "$COORDS_FILE"

read -r TAG W H SCALE STAMP <<<"$(head -n1 "$MANIFEST")"
[[ "$TAG" == CANVAS ]] || { echo "Manifest malformed"; head -n1 "$MANIFEST"; exit 1; }
echo "Canvas: ${W}x${H}px @ ${SCALE}\"/px | stamp=${STAMP}px"

# ---------- Overlap masks + stamp rectangles ----------
MASKS_FILE="masks.tsv"; : > "$MASKS_FILE"
STAMPS_FILE="stamps.tsv"; : > "$STAMPS_FILE"
awk -v STAMP="$STAMP_PX" '
NR==1{ next }
{
  n++; name[n]=$1; px[n]=$2; py[n]=$3;
  sx0[n]=px[n]-STAMP/2; sy0[n]=py[n]-STAMP/2;
  sx1[n]=sx0[n]+STAMP; sy1[n]=sy0[n]+STAMP;
}
END{
  for(i=1;i<=n;i++){
    printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", name[i], sx0[i], sy0[i], sx1[i], sy1[i], px[i], py[i] > "'$STAMPS_FILE'";
  }
  for(i=1;i<=n;i++){
    for(j=i+1;j<=n;j++){
      ox0=(sx0[i]>sx0[j]?sx0[i]:sx0[j]);
      oy0=(sy0[i]>sy0[j]?sy0[i]:sy0[j]);
      ox1=(sx1[i]<sx1[j]?sx1[i]:sx1[j]);
      oy1=(sy1[i]<sy1[j]?sy1[i]:sy1[j]);
      if(ox1>ox0 && oy1>oy0){
        rx0=ox0 - sx0[i]; ry0=oy0 - sy0[i];
        rx1=ox1 - sx0[i]; ry1=oy1 - sy0[i];
        printf "%s\t%d\t%d\t%d\t%d\n", name[i], rx0, ry0, rx1, ry1 >> "'$MASKS_FILE'";
        printf "Overlap: %s under %s\n", name[i], name[j] > "/dev/stderr";
      }
    }
  }
}' "$MANIFEST"

# ---------- Label placement (no overlap with stamps or other labels) ----------
LABELS_FILE="labels.tsv"; : > "$LABELS_FILE"
awk -v W="$W" -v H="$H" \
    -v PT="$LABEL_PT" -v CW="$LABEL_CHAR_W" -v PADX="$LABEL_PAD_X" -v PADY="$LABEL_PAD_Y" -v FIXH="$LABEL_H" \
    -v margins="$LABEL_MARGIN_LIST" '
function split_list(s, arr,   n,i,t){ n=split(s,t,/,/); for(i=1;i<=n;i++) arr[i]=t[i]+0; return n }
function ok(x0,y0,x1,y1){ return (x0>=2 && y0>=2 && x1<=W-2 && y1<=H-2) }
function inter(ax0,ay0,ax1,ay1, bx0,by0,bx1,by1){ return (ax1>bx0 && ax0<bx1 && ay1>by0 && ay0<by1) }

BEGIN{
  # read stamps
  while((getline line < "'$STAMPS_FILE'")>0){
    split(line,f,"\t"); n++; name[n]=f[1]; sx0[n]=f[2]; sy0[n]=f[3]; sx1[n]=f[4]; sy1[n]=f[5]; cx[n]=f[6]; cy[n]=f[7];
  }
  close("'$STAMPS_FILE'");
  mcount=split_list(margins, M);
  placed=0;

  for(i=1;i<=n;i++){
    txt=name[i];
    lw = int(CW*length(txt) + 2*PADX + 0.5);
    lh = FIXH;

    found=0;
    for(mi=1; mi<=mcount && !found; mi++){
      m=M[mi];
      for(pos=1; pos<=8 && !found; pos++){
        if(pos==1){ lx=int(cx[i]-lw/2); ly=int(sy0[i]-m-lh); }          # N
        else if(pos==2){ lx=int(cx[i]-lw/2); ly=int(sy1[i]+m); }        # S
        else if(pos==3){ lx=int(sx1[i]+m);    ly=int(cy[i]-lh/2); }     # E
        else if(pos==4){ lx=int(sx0[i]-m-lw); ly=int(cy[i]-lh/2); }     # W
        else if(pos==5){ lx=int(sx1[i]+m);    ly=int(sy0[i]-m-lh); }    # NE
        else if(pos==6){ lx=int(sx0[i]-m-lw); ly=int(sy0[i]-m-lh); }    # NW
        else if(pos==7){ lx=int(sx1[i]+m);    ly=int(sy1[i]+m); }       # SE
        else {           lx=int(sx0[i]-m-lw); ly=int(sy1[i]+m); }       # SW

        rx0=lx; ry0=ly; rx1=lx+lw; ry1=ly+lh;
        if(!ok(rx0,ry0,rx1,ry1)) continue;

        # must not intersect any stamp
        bad=0;
        for(s=1;s<=n;s++){ if(inter(rx0,ry0,rx1,ry1, sx0[s],sy0[s],sx1[s],sy1[s])){ bad=1; break; } }
        if(bad) continue;

        # must not intersect earlier labels
        for(k=1;k<=placed;k++){ if(inter(rx0,ry0,rx1,ry1, Lx0[k],Ly0[k],Lx1[k],Ly1[k])){ bad=1; break; } }
        if(bad) continue;

        placed++; Lname[placed]=txt; Lx0[placed]=rx0; Ly0[placed]=ry0; Lx1[placed]=rx1; Ly1[placed]=ry1; Lw[placed]=lw; Lh[placed]=lh;
        found=1;
      }
    }
    if(!found){
      # Last resort: drop farther below
      m=1200; lx=int(cx[i]-lw/2); ly=int(sy1[i]+m);
      rx0=lx; ry0=ly; rx1=lx+lw; ry1=ly+lh;
      if(ok(rx0,ry0,rx1,ry1)){ placed++; Lname[placed]=txt; Lx0[placed]=rx0; Ly0[placed]=ry0; Lx1[placed]=rx1; Ly1[placed]=ry1; Lw[placed]=lw; Lh[placed]=lh; }
    }
  }

  for(i=1;i<=placed;i++){ printf "%s\t%d\t%d\t%d\t%d\n", Lname[i], Lx0[i], Ly0[i], Lw[i], Lh[i]; }
}' > "$LABELS_FILE"

# ---------- Helper: find the cutout ----------
find_cutout() {
  local base="$1"
  for suf in "_Legacy.jpg" "_DESI_G4.jpg" ".jpg"; do
    if [[ -f "${base}${suf}" ]]; then printf "%s" "${base}${suf}"; return 0; fi
  done
  return 1
}

# ---------- Compose stamps (with overlap masking) ----------
cmd=( "$IM" "${IM_LIMITS[@]}" -size "${W}x${H}" xc:black )

mask_lines_for() { awk -v k="$1" '$1==k {print $2, $3, $4, $5}' "$MASKS_FILE"; }

while IFS=$'\t' read -r NAME PX PY; do
  [[ -z "${NAME:-}" ]] && continue
  SRC="$(find_cutout "$NAME" || true)"
  if [[ -z "$SRC" ]]; then echo "Warning: missing cutout for $NAME — skipping" >&2; continue; fi

  X=$(( PX - STAMP_PX/2 ))
  Y=$(( PY - STAMP_PX/2 ))

  MASK_LINES="$(mask_lines_for "$NAME" || true)"
  if [[ -n "${MASK_LINES}" ]]; then
    cmd+=( "(" "$SRC" -alpha on "(" -size "${STAMP_PX}x${STAMP_PX}" xc:white -fill black )
    while IFS=$' \t' read -r rx0 ry0 rx1 ry1; do
      [[ -z "${rx0:-}" ]] && continue
      (( rx0<0 )) && rx0=0; (( ry0<0 )) && ry0=0
      (( rx1>STAMP_PX )) && rx1=$STAMP_PX; (( ry1>STAMP_PX )) && ry1=$STAMP_PX
      if (( rx1>rx0 && ry1>ry0 )); then
        cmd+=( -draw "rectangle ${rx0},${ry0} ${rx1},${ry1}" )
      fi
    done <<< "${MASK_LINES}"
    cmd+=( ")" -compose CopyOpacity -composite ")" )
  else
    cmd+=( "(" "$SRC" ")" )
  fi
  cmd+=( -geometry "+${X}+${Y}" -compose over -composite )
done < <(tail -n +2 "$MANIFEST")

# ---------- Draw BIG labels onto background (PNG, crisp) ----------
while IFS=$'\t' read -r NAME LX LY LW LH; do
  [[ -z "${NAME:-}" ]] && continue
  RX=$(( LX + LW ))
  BY=$(( LY + LH ))
  TX=$(( LX + LABEL_PAD_X ))
  TY=$(( LY + LABEL_PAD_Y + LABEL_PT - 8 ))  # baseline tweak for large type
  # opaque-ish box + thick outline text for readability
  cmd+=( -fill "#000000CC" -stroke none -draw "rectangle ${LX},${LY} ${RX},${BY}" )
  cmd+=( -fill yellow -stroke black -strokewidth 10 -pointsize "${LABEL_PT}" \
         -gravity NorthWest -annotate "+${TX}+${TY}" "${NAME}" )
done < "$LABELS_FILE"

# ---------- Output (PNG, lossless) ----------
cmd+=( -colorspace sRGB "$OUTFILE" )

echo "Rendering ${OUTFILE} …"
${cmd[@]}
echo "Done → ${OUTFILE}"
