import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_violin_distances(tyr_csv_files, met_csv_files):
    """
    Plots violin plots for the distances from Trp-indole N to Tyr-O and Met-S residues across multiple families.
    
    Args:
        tyr_csv_files (list): List of CSV files for Tyr distances.
        met_csv_files (list): List of CSV files for Met distances.
    """
    # Prepare to combine Tyr distances from all datasets
    combined_tyr_df = pd.DataFrame()
    
    for file in tyr_csv_files:
        df = pd.read_csv(file)
        family_name = file.split('_')[-1].replace('.csv', '')  # Extract family name from the file name
        df['Dataset'] = family_name  # Label for the dataset
        combined_tyr_df = pd.concat([combined_tyr_df, df], ignore_index=True)

    # Prepare to combine Met distances from all datasets
    combined_met_df = pd.DataFrame()
    
    for file in met_csv_files:
        df = pd.read_csv(file)
        family_name = file.split('_')[-1].replace('.csv', '')  # Extract family name from the file name
        df['Dataset'] = family_name  # Label for the dataset
        combined_met_df = pd.concat([combined_met_df, df], ignore_index=True)

    # Create a figure with two subplots (two columns)
    fig, axs = plt.subplots(1, 2, figsize=(14, 6))

    # Violin plot for Tyr distances
    sns.violinplot(x='Dataset', y='Distance (Å)', data=combined_tyr_df, ax=axs[0], palette='muted')
    axs[0].set_title('Violin Plot of Distances from Trp-indole N to Tyr-O Residues')
    axs[0].set_xlabel('Protein Families')
    axs[0].set_ylabel('Distance (Å)')

    # Violin plot for Met distances
    sns.violinplot(x='Dataset', y='Distance (Å)', data=combined_met_df, ax=axs[1], palette='muted')
    axs[1].set_title('Violin Plot of Distances from Trp-indole N to Met-S Residues')
    axs[1].set_xlabel('Protein Families')
    axs[1].set_ylabel('Distance (Å)')

    # Adjust layout to make room for labels
    plt.tight_layout()
    plt.show()

# Example usage:
tyr_csv_files = ['output_trp_to_tyr_rho.csv', 
                 'output_trp_to_tyr_rac.csv', 
                 'output_trp_to_tyr_cdc.csv']  # Add as many Tyr CSVs as needed

met_csv_files = ['output_trp_to_met_rho.csv', 
                 'output_trp_to_met_rac.csv', 
                 'output_trp_to_met_cdc.csv']  # Add as many Met CSVs as needed

plot_violin_distances(tyr_csv_files, met_csv_files)

