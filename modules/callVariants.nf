#!/usr/bin/env nextflow

// Begin variant calling 
process callVariants {
    /*
    // Set this for local run
    container "${params.container__bcftools}"
    */

    // Set this for cluster run
    
    clusterOptions '-l select=1:ncpus=16:mem=12GB -l walltime=4:00:00 -P 12003580 -q normal'
    maxForks 20
    //publishDir "${params.output_dir}/raw_vcfs/${motifID}/", mode: 'copy'
    

    input:
    tuple val(datasetID), path(bamList)
    val(bedFiles)

    output:
    stdout
    //path('*_qualgt10.var.flt.VAF.allTFBS_${datasetID}.vcf')

    script:
    """
    #!/usr/bin/env bash

    echo \$(pwd)
    echo "list of beds:" ${bedFiles}

    # Convert the list into a string separated by spaces
    bedString="${bedFiles.join(' ')}"

    # Split the string into an array in bash
    IFS=' ' read -r -a bedArray <<< "\${bedString}"

    # Access the elements of the array
    echo "First file: \${bedArray[0]}"
    echo "Second file: \${bedArray[1]}"

    parallel -h

    mpileup_module() {
    local BEDFILES="\$1"
    local BAMLIST="\$2"
    local DATAID="\$3"
    
    echo "Running parallel job on "\$DATAID" using all 1360 motif BED files..."
    echo "\$BEDFILES" "\$BAMLIST"
    bash /mnt/runParallelMpileup.sh "\$BEDFILES" "\$BAMLIST" "\$DATAID"
    }

    # export function to be used in parallel
    export -f mpileup_module

    bedFileArray="\${bedArray[@]}"
    bamPathList="${bamList}"
    dataID="${datasetID}"

    # Determine the number of available CPU cores
    num_cores=\$(nproc)
    echo "Number of cores: "\$num_cores
    mpileup_module "\${bedFileArray}" "\${bamPathList}" "\${dataID}"

    """
}

/* # bcftools mpileup -Ou -f "${params.genome_fasta}" -T "${bedFile}" -b "${bamList}" --annotate FORMAT/AD,FORMAT/DP | bcftools call -Ou -mv | bcftools filter -i'QUAL>10' | bcftools +fill-tags - -- -t AF,VAF > "${motifID}_qualgt10.var.flt.VAF.allTFBS_${datasetID}.vcf"
*/