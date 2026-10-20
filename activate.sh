#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

if (( $# > 1 )); then
  printf 'Usage: %s [weights.tsv]\n' "$0" >&2
  exit 2
fi

weights_file=${1:-"$script_dir/anno.weights.txt"}
if [[ ! -r "$weights_file" ]]; then
  printf 'Cannot read feature weights: %s\n' "$weights_file" >&2
  exit 1
fi

awk -F '\t' '
  NR == 1 {
    if (NF != 5 || $1 != "# index" || $2 != "group" ||
        $3 != "target" || $4 != "feature" || $5 != "weight") {
      print "Invalid feature-weight header" > "/dev/stderr"
      failed = 1
      exit 1
    }
    next
  }
  {
    if (NF != 5 || $1 !~ /^[0-9]+$/ || $2 == "" || $3 == "" ||
        $4 == "" ||
        $5 !~ /^-?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/) {
      printf "Invalid feature-weight row at line %d\n", NR > "/dev/stderr"
      failed = 1
      exit 1
    }
    if ($1 != rows || seen_index[$1]++) {
      printf "Unexpected or duplicate feature index at line %d\n", NR > "/dev/stderr"
      failed = 1
      exit 1
    }
    rows++
    group = $2
    if (!(group in seen_group)) {
      seen_group[group] = 1
      group_order[++group_count] = group
    }
    total[group]++
    if ($5 + 0 > 0.01 || $5 + 0 < -0.01) active[group]++
  }
  END {
    if (failed) exit 1
    if (NR < 2 || rows == 0) {
      print "No feature weights found" > "/dev/stderr"
      exit 1
    }
    print "Feature group\tActive weights\tTotal weights\tStatus"
    all_active = 0
    for (i = 1; i <= group_count; i++) {
      group = group_order[i]
      count = active[group] + 0
      all_active += count
      printf "%s\t%d\t%d\t%s\n", group, count, total[group],
             count ? "ACTIVE" : "INACTIVE"
    }
    printf "TOTAL\t%d\t%d\t%s\n", all_active, rows,
           all_active ? "ACTIVE" : "INACTIVE"
  }
' "$weights_file"
