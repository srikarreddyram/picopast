#!/bin/zsh
# Fetch a NASA product until every file is on disk, surviving network outages.
#
#   scripts/fetch_resilient.sh PRODUCT START END PERIOD [LOGFILE]
#   e.g. caffeinate -i scripts/fetch_resilient.sh pace_rrs 2024-03-05 2026-10-31 DAY run.log
#
# Files already on disk are skipped, so a rerun only fetches what is missing.
# The retry count resets whenever new files arrived; pauses grow 1, 2, 4, 8,
# 15 min; it gives up after 8 failures in a row with no new file (~1.5 h of outage).
cd "${0:A:h}/.." || exit 1
p=$1; s=$2; e=$3; per=$4; log=${5:-/dev/stdout}
err="data/raw/.fetch_${p}_${per}_${s}.err"
count() { ls data/raw/$p 2>/dev/null | grep -c '\.box\.nc$' }
fails=0
echo "START $p $per $s..$e $(date '+%a %T')" >> $log
while true; do
  before=$(count)
  if .venv/bin/python -u .venv/bin/picopast fetch $p $s $e --period $per > /dev/null 2> $err; then
    echo "DONE $p $per $s..$e $(date '+%a %T'), $(count) files in data/raw/$p" >> $log
    exit 0
  fi
  if (( $(count) > before )); then fails=0; else fails=$((fails + 1)); fi
  echo "RETRY $p $per $s (no-progress failures: $fails) $(date '+%T'): $(grep -v -i warn $err | tail -1 | cut -c1-150)" >> $log
  (( fails >= 8 )) && { echo "FAILED $p $per $s..$e $(date '+%a %T')" >> $log; exit 1; }
  pause=$(( fails == 0 ? 30 : 60 * 2 ** (fails - 1) )); (( pause > 900 )) && pause=900
  sleep $pause
done
