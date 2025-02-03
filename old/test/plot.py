import pandas as pd
import matplotlib.pyplot as plt

def plot_distances(tyr_csv_output_file, met_csv_output_file):
    """
    Plots the distances from the CSV files using Pandas and Matplotlib in two subplots.
    
    Args:
        tyr_csv_output_file (str): Path to the CSV file for Tyr distances.
        met_csv_output_file (str): Path to the CSV file for Met distances.
    """
    # Load Tyr distances
    tyr_df = pd.read_csv(tyr_csv_output_file)
    tyr_pairs = tyr_df["PDB ID/Residue Pair"].tolist()
    tyr_distances = tyr_df["Distance (Å)"].tolist()

    # Load Met distances
    met_df = pd.read_csv(met_csv_output_file)
    met_pairs = met_df["PDB ID/Residue Pair"].tolist()
    met_distances = met_df["Distance (Å)"].tolist()

    # Create a figure with two subplots
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

    # Plot Tyr distances
    axs[0].bar(tyr_pairs, tyr_distances, color='blue')
    axs[0].set_title('Distances from Trp-indole N to Tyr-O Residues')
    axs[0].set_xlabel('PDB ID / Residue Pair')
    axs[0].set_ylabel('Distance (Å)')
    axs[0].tick_params(axis='x', rotation=45)
    axs[0].grid(axis='y')  # Add gridlines for better readability

    # Plot Met distances
    axs[1].bar(met_pairs, met_distances, color='orange')
    axs[1].set_title('Distances from Trp-indole N to Met-S Residues')
    axs[1].set_xlabel('PDB ID / Residue Pair')
    axs[1].set_ylabel('Distance (Å)')
    axs[1].tick_params(axis='x', rotation=45)
    axs[1].grid(axis='y')  # Add gridlines for better readability

    # Adjust layout to make room for x-axis labels
    plt.tight_layout()
    plt.show()

# Example usage:
tyr_csv_output_file = 'output_trp_to_tyr.csv'  # Output CSV for Tyr distances
met_csv_output_file = 'output_trp_to_met.csv'  # Output CSV for Met distances
plot_distances(tyr_csv_output_file, met_csv_output_file)
