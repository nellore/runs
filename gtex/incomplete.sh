#!/usr/bin/env bash
# Computes stats on incompletely downloaded samples. First command-line parameter should be GTEx download dir with batch_* subdirectories.
GTEXDIR=$1
pypy incomplete.py --gtex-dir $GTEXDIR | (read -r; printf "%s\n" "$REPLY"; sort -k7,7nr) >incomplete.tsv
