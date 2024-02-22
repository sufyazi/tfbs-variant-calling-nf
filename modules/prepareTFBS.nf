#!/usr/bin/env nextflow

process prepareTFBS {
    // Run inside a container with Python/Pandas installed
    container "${params.container__bcftools}"

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