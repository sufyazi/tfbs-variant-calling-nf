#!/usr/bin/env nextflow

// Begin variant calling 
process callVariants {
    
    cpus 8

    // Set this for local run
    container "${params.container__bcftools}"

    // Set this for cluster run
    /*
    clusterOptions '-l select=1:ncpus=64:mem=120GB -l walltime=4:00:00 -P 12003580 -q normal'
    maxForks 20
    publishDir "${params.output_dir}/raw_vcfs/${motifID}/", mode: 'copy'
    */

    input:
    tuple val(datasetID), path(bamList)
    each bedFile

    output:
    stdout
    //path "${motifID}_qualgt10.var.flt.VAF.allTFBS_${datasetID}.vcf"

    script:
    """
    FILE=${bedFile}
    motifID=\$(basename \${FILE})
    echo hello world, ${datasetID}, \${motifID%%_TOBIAS_TF_binding_sites-sorted.bed}, ${bamList}

    """
}

/* # bcftools mpileup -Ou -f "${params.genome_fasta}" -T "${bedFile}" -b "${bamList}" --annotate FORMAT/AD,FORMAT/DP | bcftools call -Ou -mv | bcftools filter -i'QUAL>10' | bcftools +fill-tags - -- -t AF,VAF > "${motifID}_qualgt10.var.flt.VAF.allTFBS_${datasetID}.vcf"
*/