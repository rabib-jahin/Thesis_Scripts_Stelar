import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Assuming your data is in a CSV file named '../Dataset/final.csv'
data = pd.read_csv('../Dataset/final.csv')

# Splitting data by 'folder'
folders = data['folder'].unique()

# Creating a bar chart for each folder
for folder in folders:
    folder_data = data[data['folder'] == folder]

    # Getting unique inner_folders for the x-axis
    inner_folders = folder_data['inner_folder'].unique()

    # Preparing the data for plotting
    summary_methods = folder_data['summary_method'].unique()
    rooting_methods = folder_data['rooting_method'].fillna('None').unique()

    # Creating a bar for each summary_method-rooting_method combination
    bar_width = 0.2
    x = np.arange(len(inner_folders))

    fig, ax = plt.subplots()
    for i, summary_method in enumerate(summary_methods):
        for j, rooting_method in enumerate(rooting_methods):
            # Initializing y values for all inner_folders
            y = np.zeros(len(inner_folders))
            for k, inner_folder in enumerate(inner_folders):
                # Get average_rf for specific combination, default to 0 if missing
                y[k] = folder_data[(folder_data['summary_method'] == summary_method) & 
                                   (folder_data['rooting_method'].fillna('None') == rooting_method) & 
                                   (folder_data['inner_folder'] == inner_folder)]['average_rf'].fillna(0).mean()
            
            ax.bar(x + (i*len(rooting_methods)+j) * bar_width, y, width=bar_width, label=f'{summary_method}-{rooting_method}')

    # Adding labels and title
    ax.set_xlabel('Inner Folder')
    ax.set_ylabel('Average RF')
    ax.set_title(f'Average RF for {folder}')
    ax.set_xticks(x + bar_width * (len(summary_methods) * len(rooting_methods)) / 2)
    ax.set_xticklabels(inner_folders)
    ax.legend()

    # Display the plot
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()
