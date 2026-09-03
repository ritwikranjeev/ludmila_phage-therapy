#!/usr/bin/env bash
# Supplementary Figure 1 — clinical timeline, first presentation through
# post-treatment follow-up. One script, one composite output, no panels.
#
#   ./run.sh
#
# Time is in days relative to the first intravenous vancomycin dose (day 0),
# the same anchor as Figure 3. The input carries no calendar dates.
#
# Entries with no recorded end date are drawn 14 days wide. That width is an
# assumption; the script lists every affected entry in logs/supp01.log. To make
# the assumption visible in the figure instead, shrink it:
#
#   S1_SPAN_DAYS=0 ./run.sh
set -uo pipefail
cd "$(dirname "$0")"
mkdir -p logs

echo "Supplementary Figure 1"
if FIG_DATA="$PWD/data" FIG_OUT="$PWD/output" \
   Rscript scripts/clinical_timeline.R > logs/supp01.log 2>&1; then
  printf '  OK    clinical_timeline.R -> output/\n'
else
  printf '  FAIL  clinical_timeline.R -> logs/supp01.log\n'
  exit 1
fi
