#!/usr/bin/env bash
# shellcheck disable=SC1091,SC2034

set -euo pipefail

echo "Finding input file paths..."

# check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: findPaths.sh <input_folder> <dataset_id_list_one_col>"
    exit 1
fi

# set up variables
bam_rootdir=$1
dataset_ids=$2

# loop through the dataset IDs
readarray -t data_ids < "$dataset_ids"

for data_id in "${data_ids[@]}"; do
    echo "Dataset ID: $data_id"
    # find all bam files in the dataset directory
    readarray -t bams < <(find "$bam_dir"/"$data_id" \( -name "*nodup.no_chrM_MT.bam" -o -name "*nodup.rep-merged.bam" \) -type f)
    echo "Number of bam files: " "${#bams[@]}"
    # check if bam file variable is not empty
    if [ -n "${bams[*]}" ]; then
        echo "Found bam files. Proceeding..."
        echo "${bams[@]}"

        echo "Submitting PBS jobs sequentially..."
        # # submit job
        # qsub -v TF_LIST="$3",BAM_INP="${bam_files}",OUT_DIR="${outfile}",ID="${data_id}" /home/users/ntu/suffiazi/scripts/gatk-workflow-scripts/scripts/customs/run-bcftools_mpileup_parallel-v2.pbs
    else
        echo "No bam files found in $bam_dir/$data_id"
        exit 1
    fi
done


echo "DONE!"