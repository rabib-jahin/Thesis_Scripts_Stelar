#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Path to the CSV file
csv_file_path = 'rf_average.csv'  # Replace with your actual CSV file path

# Read the data from the CSV file into a Pandas DataFrame
df = pd.read_csv('rf_average.csv')

# Assume the CSV has columns 'folder', 'inner_folder', 'method', 'average_rf', 'std_rf'
# which correspond to the folder names, method names, average RF scores, and standard deviations, respectively.

# Plot settings
bar_width = 0.2  # Width of the bars
opacity = 0.8    # Opacity of the bars

# Set up the figure and axes for the bar chart
fig, ax = plt.subplots(figsize=(10, 6))

# Generate a color for each method
colors = [plt.get_cmap('tab10')(i) for i in range(len(df['method'].unique()))]

# Group the data by 'inner_folder' to create grouped bar chart
for i, (name, group) in enumerate(df.groupby('inner_folder')):
    # Calculate bar positions
    index = np.arange(len(group)) * len(df['method'].unique()) + i
    # Plot the bars
    bars = ax.bar(index, group['average_rf'], bar_width,
                  alpha=opacity, color=colors[i % len(colors)],
                  label=name, yerr=group['std_rf'],
                  error_kw={'ecolor': 'black', 'capsize': 5})

# Final plot adjustments
ax.set_xlabel('Method')
ax.set_ylabel('Average RF Score')
ax.set_title('Average RF Scores by Method with Error Bars')
ax.set_xticks(np.arange(len(df['method'].unique())) * len(df['method'].unique()) + bar_width)
ax.set_xticklabels(df['method'].unique())
ax.legend(title='Inner Folder')

# Display the plot
plt.tight_layout()
plt.show()
