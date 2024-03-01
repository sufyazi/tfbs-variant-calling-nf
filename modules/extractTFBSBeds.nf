#!/usr/bin/env nextflow

// Extract TFBS from TOBIAS footprint motif matrices
process extractTFBSBeds {

    // Set this for local run
    container "${params.container__bcftools}"

    // Set this for cluster run
    /*
    clusterOptions '-l select=1:ncpus=1:mem=16GB -l walltime=4:00:00 -P 12003580 -q normal'
    maxForks 20
    publishDir "${params.output_dir}/sorted_beds/", mode: 'copy', overwrite: true
    */
    
    input:
        tuple val(motifid), path(matrix)

    output:
        path("${motifid}_TOBIAS_TF_binding_sites-sorted.bed"), emit: sorted_bed

    script:
    """
    extractTFBSRegions.py $matrix $motifid && sortBeds.sh "${motifid}_TOBIAS_TF_binding_sites-unsorted.bed"
    """
}