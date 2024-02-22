#!/usr/bin/env nextflow

// All of the default parameters are being set in `nextflow.config`
// Import sub-workflows
include { prepareTFBS } from './modules/prepareTFBS'


// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run main.nf <ARGUMENTS>

Required Arguments:

  Input Data:
  --fpscore_matrix        Folder containing directories of datasets containing BAM files ending with .bam, alongside a corresponding .bam.bai file

  TFBS Matrices:
  --tfbs_prefix_list        Motif ID prefix for the TFBS matrices

  Output Location:
  --output_dir       Folder for output files; must contain the subfolders `sorted_beds` and `raw_vcfs`
    """.stripIndent()
}


// Main workflow
workflow {

    // Show help message if the user specifies the --help flag at runtime
    // or if any required params are not provided
    if ( params.help || params.output_dir == false || params.genome_fasta == false ){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    }

    // // The user should specify --bam_folder AND --dataset_id_list
    // if ( ! params.bam_folder || ! params.dataset_id_list ){
    //     log.info"""
    //     User must specify --bam_folder AND --dataset_id_list
    //     """.stripIndent()
    //     // Exit out and do not run anything else
    //     exit 1
    // }

    // if ( params.bam_folder && params.dataset_id_list){
    //     log.info"""
    //     User has specified --bam_folder AND --dataset_id_list. Proceeding with workflow...
    //     """.stripIndent()
        
    // Preparing the input data
    log.info "Preparing the input data..."

    // Set up a channel to grab all the matrix files in the input folder
    in_matrix_ch = Channel.fromPath(params.fpscore_matrix).view()
    
    // Extract the prefix from the input files
    motif_matrices_ch = in_matrix_ch.map{ file -> [file, file.baseName.replaceAll("_tfbs_merged_matrix-full", "")] }.view()

    // Prepare the TFBS matrices
    prepareTFBS(motif_matrices_ch)


        // create a channel of dataset IDs from the dataset_id_list file
        // dataset_id_ch  = Channel
        //                         .fromPath(params.dataset_id_list)
        //                         .splitText()
        //                         .map{ it.trim() }

        // // Define the pattern which will be used to find the FASTQ files
        // fastq_pattern = "${params.fastq_folder}/*_R{1,2}*fastq.gz"

        // // Set up a channel from the pairs of files found with that pattern
        // fastq_ch = Channel
        //     .fromFilePairs(fastq_pattern)
        //     .ifEmpty { error "No files found matching the pattern ${fastq_pattern}" }
        //     .map{
        //         [it[0], it[1][0], it[1][1]]
        //     }


    // // Perform quality trimming on the input 
    // quality_wf(
    //     fastq_ch
    // )
    // // output:
    // //   reads:
    // //     tuple val(specimen), path(read_1), path(read_2)

    // // Align the quality-trimmed reads to the reference genome
    // align_wf(
    //     quality_wf.out.reads,
    //     file(params.genome_fasta)
    // )
    // // output:
    // //   bam:
    // //     tuple val(specimen), path(bam)

}