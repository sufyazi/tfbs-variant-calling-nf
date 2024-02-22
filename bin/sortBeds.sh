#!/usr/bin/env bash

set -euo pipefail

bed_file=$1

echo "Sorting bed file"
sort -k1,1V -k2,2n "$bed_file" > "${bed_file%-unsorted.bed}-sorted.bed"

echo DONE!
