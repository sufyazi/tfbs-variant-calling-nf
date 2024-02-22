#!/usr/bin/env python3

import os
import sys
import pandas as pd
import pyarrow.dataset as ds

# get arguments
in_matrix = sys.argv[1]
tf_prefix = sys.argv[2]

# Create a dataset object from the parquet file
dataset = ds.dataset(in_matrix, format="parquet")
# Get the list of columns in the dataset
column_names = dataset.schema.names
# Select columns by index
selected_columns = [column_names[i] for i in [0, 1, 2, 3]]

df = dataset.to_table(columns=selected_columns).to_pandas()

# construct the output file path
output = os.path.join(tf_prefix + "_TOBIAS_TF_binding_sites-unsorted.bed")

# remove header for bedtools and then save to file
df.to_csv(output, sep="\t", header=False, index=False)
