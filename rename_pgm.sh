#! /usr/bin/env bash
# exit on fail
# set -e

cd workloadc

for f in *.txt; do mv "$f" "$(echo "$f" | sed s/PGM_4_/PGM_9_/)"; done
for f in *.txt; do mv "$f" "$(echo "$f" | sed s/PGM_16_/PGM_7_/)"; done
for f in *.txt; do mv "$f" "$(echo "$f" | sed s/PGM_32_/PGM_6_/)"; done
for f in *.txt; do mv "$f" "$(echo "$f" | sed s/PGM_64_/PGM_5_/)"; done
for f in *.txt; do mv "$f" "$(echo "$f" | sed s/PGM_256_/PGM_4_/)"; done
for f in *.txt; do mv "$f" "$(echo "$f" | sed s/PGM_1024_/PGM_3_/)"; done
for f in *.txt; do mv "$f" "$(echo "$f" | sed s/PGM_2048_/PGM_2_/)"; done
for f in *.txt; do mv "$f" "$(echo "$f" | sed s/PGM_4096_/PGM_1_/)"; done
for f in *.txt; do mv "$f" "$(echo "$f" | sed s/PGM_8192_/PGM_0_/)"; done
