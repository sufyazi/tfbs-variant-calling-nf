#!/usr/bin/env nextflow

// All of the default parameters are being set in `nextflow.config`
// Import sub-workflows
include { extractTFBSBeds } from './modules/extractTFBSBeds'
include { generateBAMPaths } from './modules/generateBAMPaths'

// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run main.nf <--OPTION NAME> <ARGUMENT>

Required Arguments:
    -profile                      Nextflow profile to use for the subworkflows. Specify `bedfile_extraction` to run the TFBS bed file generation subworkflow. Omit this Nextflow option if default profile (run main workflow only) is desired. 
    --run_mode                    Set the run mode for the workflow. Options are `all` for all available dataset IDs in the base directory [--bam_dir] or `subset` for just the IDs specified in the list provided with [--dataset_id_list]. Default is `subset`
    --help                        Print this help message and exit

  Input paths:
    --bam_dir                     Path to base directory where directories of datasets containing BAM files and their respective bai index files are located [MANDATORY]
    --genome_fasta                Path to the reference genome fasta file [MANDATORY]
    --fpscore_matrix              Path to base directory where the merged footprint score matrix files of all studied motifs are located. Required if `-profile bedfile_extraction` is set

  Input manifests:
    --tfbs_prefix_list            Path to a list of motif ID prefixes for the TFBS matrices. Required if [-profile bedfile_extraction] is set
    --dataset_id_list             Path to a list of dataset IDs to work on. [MANDATORY as --run_mode is set to `subset` by default]

  Output path:
    --output_dir                  Directory path for output VCF files [MANDATORY]
    """.stripIndent()
}


// Main workflow
workflow {
    // Show help message if the user specifies the --help flag at runtime
    // or if any required params are not provided
    if ( params.help || params.bam_dir == false || params.output_dir == false || params.genome_fasta == false ){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    } else {
        log.info "Proceeding with workflow..."

        // Preparing the input data
        log.info "Preparing the input data..."

        // Check if the bedfile_extraction profile is set
        if ( params.subworkflow ){
            log.info "The <bedfile_extraction> profile has been set. Running the TFBS bed file generation subworkflow..."
            // Check if fpscore_matrix is provided
            if (params.fpscore_matrix == false){
                log.error "The [--fpscore_matrix] parameter is required to run TFBS bed file generation."
                exit 1
            }
            // Set up a channel to grab all the matrix files in the input folder
            in_matrix = Channel.fromPath(params.fpscore_matrix)//.view()
            // Extract the prefix from the input files and return a tuple of the file and the prefix
            motifMatrix_ch = in_matrix.map{ file -> [file, file.baseName.replaceAll("_tfbs_merged_matrix-full", "")] }//.view()
            // Also create a channel just of the prefixes
            motifPrefix_ch = motifMatrix_ch.map{ it[1] }//.view()
        
            // Run the sub-workflow to extract the TFBS as sorted bed files
            // Extract the TF footprint regions (TFBS) from the input fps matrix files
            extractTFBSBeds(motifMatrix_ch)

        } else {
            log.info "Skipping the extraction of TFBS as bed files..."

            // Check if tfbs_prefix_list is provided
            if (params.tfbs_prefix_list == false){
                log.error "The [--tfbs_prefix_list] parameter is required to run the main workflow if TFBS bed file generation is skipped."
                exit 1
            } else {
                // Set up a channel to grab all the motif ID prefixes provided in the input list
                motifPrefix_ch = Channel.fromPath(params.tfbs_prefix_list).splitText().map { it.trim() }//.view()
            }
        }
        
        // Check if run_mode is set
        if (params.run_mode == "subset"){
            log.info "The <subset> parameter for [--run_mode] has been set. Checking if input dataset ID list is provided..."
            // Check if dataset_id_list is provided
            if (params.dataset_id_list == false){
                log.error "The [--dataset_id_list] option is required to run the workflow if [--run_mode] is set to 'subset'."
                exit 1
            } else {
                log.info "A dataset ID list has been provided. Extracting BAM directory paths only for the specified IDs..."
                // Extract bam directory paths
                datasetIDs = Channel.fromPath(params.dataset_id_list).splitText().map { it.trim() }
                datasetIDPaths_ch = datasetIDs.map { id -> [id, file("${params.bam_dir}/${id}")] }//.view()
            }
        }
        if (params.run_mode == "all"){
            log.info "The <all> parameter for [--run_mode] has been set. Extracting all available dataset IDs within the input BAM directory..."
            // Extract all unique dataset IDs in the input bam folder and the path to the dataset ID
            datasetIDPaths_ch = Channel.fromPath("${params.bam_dir}/*", type: 'dir').map { dir -> [dir.name, dir] }//.view()
        }

        // Set up a channel to grab all the bam files for each dataset ID
        datasetIDBams_ch = datasetIDPaths_ch.map { id, path -> 
                                            def bams = files("${path}/**.bam")
                                            return [id, bams]
                                            }//.view()

        // Generate a list of all the bam files for all the dataset IDs 
        generateBAMPaths(datasetIDBams_ch)
        
        // now we can run the variant-calling process
        /* variantCallmpileup(
            generateBAMPaths.out.bam
        )*/

    }
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}