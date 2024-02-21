#!/usr/bin/env python3

import os
import sys
import pandas as pd
import pyarrow.parquet as pq
import pyarrow.dataset as ds

# Read arguments
input_filepath = sys.argv[1]
tf_prefix = sys.argv[2]

# Create a dataset object from the parquet file
dataset = ds.dataset(input_filepath, format="parquet")
df = pq.read_table(input_filepath).to_pandas()
            # process the matrix file; grab the first 4 cols, add # to the header line for bedtools, then use sed to sort the bed file by chr and start position, ignoring the header line (1q)
            awk -v OFS='\t' '{print $1, $2, $3, $4}' "${matrix}" | sed '1s/^/#/' | (sed -u 1q; sort -k1,1V -k2,2n) > "$output_dir"/"${prefix}"_diffmode_TCGA-BRCA_fpscore_regions.bed
        fi
    fi
done < "$input_prefix"
