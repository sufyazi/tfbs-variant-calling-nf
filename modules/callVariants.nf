#!/usr/bin/env nextflow

// Begin variant calling 
process callVariants {
    /*
    // Set this for local run
    container "${params.container__bcftools}"
    */

    // Set this for cluster run
    
    clusterOptions '-l select=1:ncpus=64:mem=120GB -l walltime=4:00:00 -P 12003580 -q normal'
    maxForks 20
    publishDir "${params.output_dir}/raw_vcfs/${datasetID}/", mode: 'copy'
    

    input:
    tuple val(datasetID), path(bamList)
    val(bedFiles)

    output:
    path('*_qualgt10.var.flt.VAF.allTFBS_*.txt')

    script:
    """
    #!/usr/bin/env bash

    # Convert the list into a string separated by spaces
    bedString="${bedFiles.join(' ')}"

    # Split the string into an array in bash
    IFS=' ' read -r -a bedArray <<< "\${bedString}"

    # Count number of elements in array
    echo "No. of BED files: " "\${#bedArray[@]}"

    mpileup_module() {
    local BEDFILES="\$1"
    local BAMLIST="\$2"
    local DATAID="\$3"
    local FASTAFILE="\$4"
    
    bash /mnt/runParallelMpileup.sh "\$BEDFILES" "\$BAMLIST" "\$DATAID" "\$FASTAFILE"
    }

    # export function to be used in parallel
    export -f mpileup_module

    # determine the number of available CPU cores
    num_cores=\$(nproc)
    echo "Number of cores: "\$num_cores

    # set arguments
    bedFileArray="\${bedArray[@]}"
    bamPathList="${bamList}"
    dataID="${datasetID}"

    # Use GNU Parallel to run the commands in parallel using all available cores
    parallel -j "\$num_cores" mpileup_module ::: "\${bedArray[@]}" ::: "\${bamPathList}" ::: "\${dataID}" ::: "${params.genome_fasta}"

    """
}