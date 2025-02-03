import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_violin_distances(tyr_csv_files, met_csv_files):
    """
    Plots violin plots for the distances from Trp-indole N to Tyr-O and Met-S residues across multiple families,
    split by TRP residue numbers (excluding PDB codes).
    
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
        
        # Combine dataset name and TRP residue number
        df['Residue ID'] = df.apply(lambda x: f"{family_name}-{x['TRP Residue']}", axis=1)
        combined_tyr_df = pd.concat([combined_tyr_df, df], ignore_index=True)

    # Prepare to combine Met distances from all datasets
    combined_met_df = pd.DataFrame()
    
    for file in met_csv_files:
        df = pd.read_csv(file)
        family_name = file.split('_')[-1].replace('.csv', '')  # Extract family name from the file name
        df['Dataset'] = family_name  # Label for the dataset
        
        # Combine dataset name and TRP residue number
        df['Residue ID'] = df.apply(lambda x: f"{family_name}-{x['TRP Residue']}", axis=1)
        combined_met_df = pd.concat([combined_met_df, df], ignore_index=True)

    # Create a figure with two subplots (two columns)
    fig, axs = plt.subplots(1, 2, figsize=(16, 6), sharey=True)

    # Violin plot for Tyr distances
    sns.violinplot(x='Residue ID', y='Distance (Å)', data=combined_tyr_df, ax=axs[0], palette='muted')
    axs[0].set_title('Violin Plot of Distances from Trp-indole N to Tyr-O Residues')
    axs[0].set_xlabel('TRP Residue Numbers by Family')
    axs[0].set_ylabel('Distance (Å)')
    axs[0].tick_params(axis='x', rotation=45)

    # Violin plot for Met distances
    sns.violinplot(x='Residue ID', y='Distance (Å)', data=combined_met_df, ax=axs[1], palette='muted')
    axs[1].set_title('Violin Plot of Distances from Trp-indole N to Met-S Residues')
    axs[1].set_xlabel('TRP Residue Numbers by Family')
    axs[1].set_ylabel('Distance (Å)')
    axs[1].tick_params(axis='x', rotation=45)

    # Adjust layout to make room for labels
    plt.tight_layout()
    #plt.savefig("dist_vio.pdf")
    plt.show()

# Example usage:
# Tyr CSVs for rho, rac, cdc
tyr_csv_files = ['output_trp_to_tyr_rho.csv', 
                 'output_trp_to_tyr_rac.csv', 
                 'output_trp_to_tyr_cdc.csv']

# Met CSVs for rho, rac, cdc
met_csv_files = ['output_trp_to_met_rho.csv', 
                 'output_trp_to_met_rac.csv', 
                 'output_trp_to_met_cdc.csv']

plot_violin_distances(tyr_csv_files, met_csv_files)
