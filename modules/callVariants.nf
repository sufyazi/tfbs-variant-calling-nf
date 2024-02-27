#!/usr/bin/env nextflow

// Begin variant calling 
process callVariants {
    
    clusterOptions '-l select=1:ncpus=1:mem=400GB -l walltime=18:00:00 -P 12003580 -q normal'
    
    //publishDir "${params.output_dir}/raw_vcfs/", mode: 'move'

    input:
    tuple val(motifID), path(bedFile), val(datasetID), path(bamList)

    output:
    path "${motifID}_qualgt10.var.flt.VAF.allTFBS_${datasetID}.vcf"

    script:
    """
    echo "Starting BCFTools from within the container for ${motifID} on dataset ${datasetID}..."
    
    bcftools mpileup -Ou -f "${params.genome_fasta}" -T "${bedFile}" -b "${bamList}" --annotate FORMAT/AD,FORMAT/DP | bcftools call -Ou -mv | bcftools filter -i'QUAL>10' | bcftools +fill-tags - -- -t AF,VAF > "${motifID}_qualgt10.var.flt.VAF.allTFBS_${datasetID}.vcf"
    
    """
}

/* 
*/