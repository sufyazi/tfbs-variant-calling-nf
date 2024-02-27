#!/bin/bash

timestamp=$(date +%Y%m%d-%H.%M)

# Run the workflow on the test data, and write the output to output/
nextflow \
    run \
    main.nf \
    --run_mode "subset" \
    --dataset_id_list "input-dir-cluster/manifests/dataset-manifest-test.txt" \
    --fpscore_matrix "input-dir-cluster/test-matrix" \
    --genome_fasta "input-dir-cluster/refs/GRCh38_no_alt_GCA_000001405.15.fasta.gz" \
    --output_dir "output-dir-cluster" \
    -c nextflow_cluster.config \
    -with-dag workflow-reports/dag-${timestamp}.html \
    -with-report workflow-reports/report-${timestamp}.html \
    -with-timeline workflow-reports/timeline-${timestamp}.html
    # -resume \
