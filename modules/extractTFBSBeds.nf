#!/usr/bin/env nextflow

// Extract TFBS from TOBIAS footprint motif matrices
process extractTFBSBeds {

    // container "${params.container__bcftools}"
    clusterOptions '-l select=1:ncpus=1:mem=10GB -l walltime=4:00:00 -P 12003580 -q normal'

    publishDir "${params.output_dir}/sorted_beds/", mode: 'copy', overwrite: true

    input:
        tuple path(matrix), val(motifid)

    output:
        path("${motifid}_TOBIAS_TF_binding_sites-sorted.bed")

    script:
    """
    extractTFBSRegions.py $matrix $motifid && sortBeds.sh "${motifid}_TOBIAS_TF_binding_sites-unsorted.bed"
    """
}