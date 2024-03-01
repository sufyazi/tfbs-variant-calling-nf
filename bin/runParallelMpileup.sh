#!/usr/bin/env bash

set -euo pipefail

bedfile=$1
bamlist=$2
datasetid=$3
fasta=$4

basefile=$(basename "${bedfile}")

echo "Parallel job to execute with mpileup: $bedfile" "$bamlist" "$datasetid" "$fasta"

bcftools mpileup -Ou -f "${fasta}" -T "${bedfile}" -b "${bamlist}" --annotate FORMAT/AD,FORMAT/DP | bcftools call -Ou -mv | bcftools filter -i'QUAL>10' | bcftools +fill-tags - -- -t AF,VAF > "${basefile%%_TOBIAS_TF_binding_sites-sorted.bed}_qualgt10.var.flt.VAF.allTFBS_${datasetid}.vcf"

