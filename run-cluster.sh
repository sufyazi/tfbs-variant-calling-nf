#!/bin/bash

# Run the workflow on the test data, and write the output to output/
nextflow \
    run \
    -resume \
    main.nf \
    --run_mode "subset" \
    --dataset_id_list "test_data/dataset-manifest-heme-cancers.txt" \
    --genome_fasta "temp_dir-cluster/refs/GRCh38_no_alt_GCA_000001405.15.fasta.gz" \
    --output_dir "test_data/output" \
    -c nextflow_cluster.config \
    -with-dag
    # 
