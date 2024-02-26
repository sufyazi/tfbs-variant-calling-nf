#!/bin/bash

# Run the workflow on the test data, and write the output to output/
nextflow \
    run \
    main.nf \
    -profile bedfile_extraction \
    --run_mode "all" \
    --bam_dir "test_data/bam_dir" \
    --genome_fasta "temp_dir/GRCh38_no_alt_GCA_000001405.15.fasta.gz" \
    --output_dir "test_data/output" # \
    # 
    # 
    # -with-report \
    # -resume
