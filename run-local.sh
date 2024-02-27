#!/bin/bash

# Run the workflow on the test data, and write the output to output/
nextflow \
    run \
    main.nf \
 
    # -with-report \
    # -resume
