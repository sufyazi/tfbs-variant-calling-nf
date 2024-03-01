#!/bin/bash

# Run the workflow on the test data, and write the output to output/
nextflow \
    -log test-local.log \
    run \
    -resume \
    main.nf \
    --output_dir "output-dir-local" \
    -c nextflow_local.config
