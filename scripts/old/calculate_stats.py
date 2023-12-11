#!/usr/bin/env python3
import sys
import numpy as np
import pandas as pd

# Read the file path from command line arguments
rf_scores_file = sys.argv[1]
rf_average_file = sys.argv[2]

# Read the RF scores into a pandas DataFrame
rf_df = pd.read_csv(rf_scores_file, header=None, names=['folder', 'inner_folder', 'method', 'rf'])

# Calculate mean and standard deviation
stats_df = rf_df.groupby(['folder', 'inner_folder', 'method']).agg(['mean', 'std']).reset_index()
stats_df.columns = ['folder', 'inner_folder', 'method', 'average_rf', 'std_rf']

# Write the averages and standard deviations to a CSV file
stats_df.to_csv(rf_average_file, index=False)

print(f"RF score statistics calculated and written to {rf_average_file}")
