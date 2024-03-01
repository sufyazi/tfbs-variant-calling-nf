#!/usr/bin/env nextflow

// Extract the BAM file paths from the input folder
process generateBAMPaths {
    /*
    // Set this for local run
    container "${params.container__bcftools}"
    */

    // Set this for cluster run
    clusterOptions '-l select=1:ncpus=1:mem=16GB -l walltime=4:00:00 -P 12003580 -q normal'
    maxForks 10
    
    
    input:
    tuple val(datasetID), val(bams)

    output:
    tuple val(datasetID), path("${datasetID}.bamlist.txt")

    script:
    def bams_formatted = bams.join("\n")
    """
    echo "${bams_formatted}" > ${datasetID}.bamlist.txt
    """
}