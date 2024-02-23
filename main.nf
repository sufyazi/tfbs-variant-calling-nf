#!/usr/bin/env nextflow

// All of the default parameters are being set in `nextflow.config`
// Import sub-workflows
include { extractTFBSBeds } from './modules/extractTFBSBeds'
include { generateBAMPaths } from './modules/generateBAMPaths'

// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run main.nf <ARGUMENTS>

Required Arguments:

  Input Data:
  --fpscore_matrix              Folder containing directories of datasets containing BAM files ending with .bam, alongside a corresponding .bam.bai file

  TFBS Matrices:
  --tfbs_prefix_list            Motif ID prefix for the TFBS matrices

  Output Location:
  --output_dir                  Folder for output files; must contain the subfolders `sorted_beds` and `raw_vcfs`
    """.stripIndent()
}


// Main workflow
workflow {

    // Show help message if the user specifies the --help flag at runtime
    // or if any required params are not provided
    if ( params.help || params.output_dir == false || params.genome_fasta == false || params.bam_dir == false ){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    } else {
        log.info "Proceeding with workflow..."

        // Preparing the input data
        log.info "Preparing the input data..."

        if ( params.extract_tfbs_as_beds ){
            // Check if fpscore_matrix is provided
            if (params.fpscore_matrix == false){
                log.error "The --fpscore_matrix parameter is required to extract TFBS as bed files"
                exit 1
            }
            // Set up a channel to grab all the matrix files in the input folder
            in_matrix = Channel.fromPath(params.fpscore_matrix)
            // Extract the prefix from the input files and return a tuple of the file and the prefix
            motifMatrix_ch = in_matrix.map{ file -> [file, file.baseName.replaceAll("_tfbs_merged_matrix-full", "")] }.view()
            // Also create a channel just of the prefixes
            motifPrefix_ch = motifMatrix_ch.map{ it[1] }
        
            // Run the sub-workflow to extract the TFBS as sorted bed files
            // Extract the TF footprint regions (TFBS) from the input fps matrix files
            extractTFBSBeds(motifMatrix_ch)
            
        } else {
            log.info "Skipping the extraction of TFBS as bed files..."

            // Check if tfbs_prefix_list is provided
            if (params.tfbs_prefix_list == false){
                log.error "The --tfbs_prefix_list parameter is required to run the main workflow if TFBS bed file generation is skipped."
                exit 1
            } else {
                // Set up a channel to grab all the matrix files in the input folder
                motifPrefix_ch = Channel.fromPath(params.tfbs_prefix_list).splitText().map { it.trim() }//.view()
            }    
        }

        // test if param.dataset_id_list is provided or not
        if (params.dataset_id_list){
            log.info "A dataset ID list has been provided. Extracting BAM directory paths only for the specified IDs..."
            // If the user has provided a list of dataset IDs, use that to extract bam directory paths
            datasetIDs = Channel.fromPath(params.dataset_id_list).splitText().map { it.trim() }
            datasetIDPaths_ch = datasetIDs.map { id -> [id, file("${params.bam_dir}/${id}")] }//.view()

        } else {
            log.info "No dataset ID list provided. Extracting all available dataset IDs within the input BAM directory..." 
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
    }

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