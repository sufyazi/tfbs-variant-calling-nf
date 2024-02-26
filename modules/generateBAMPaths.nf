#!/usr/bin/env nextflow

// Extract the BAM file paths from the input folder
process generateBAMPaths {
    
    // container "${params.container__bcftools}"
    clusterOptions '-l select=1:ncpus=1:mem=10GB -l walltime=4:00:00 -P 12003580 -q normal'
    
    input:
    tuple val(datasetID), val(bams)

    output:
    path "${datasetID}.bamlist.txt", emit: bamlist

    script:
    def bams_formatted = bams.join("\n")
    """
    echo "${bams_formatted}" > ${datasetID}.bamlist.txt
    """
}