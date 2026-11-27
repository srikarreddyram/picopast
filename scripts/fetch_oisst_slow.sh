#!/bin/zsh
# Fetch missing OISST box years one at a time; on NOAA PSL rate-limiting, wait and retry.
#   scripts/fetch_oisst_slow.sh CONFIG RAWDIR FIRST LAST LOG
cd "${0:A:h}/.." || exit 1
cfg=$1; raw=$2; first=$3; last=$4; log=$5
for y in $(seq $first $last); do
  [[ -f $raw/oisst/oisst_v2.1_box_$y.nc && $y -lt $(date +%Y) ]] && continue
  until .venv/bin/python -W ignore -m picopast.cli --config $cfg fetch-open --start-year $y --end-year $y >> $log 2>&1; do
    echo "RETRY $y $(date '+%T'): waiting 20 min" >> $log; sleep 1200
  done
  echo "DONE $y $(date '+%T')" >> $log
done
echo "ALL DONE $(date '+%a %T')" >> $log
