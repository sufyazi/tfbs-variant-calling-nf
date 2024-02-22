#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Process the bam file paths and store them in a channel
process findPaths {
    
    input:
    tuple path(bam_file), val(tf_prefix)

    output:
    path "fastqc/*.zip", emit: zip
    path "fastqc/*.html", emit: html

    script:
    template 'fastqc.sh'

}