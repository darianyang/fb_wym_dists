import MDAnalysis as mda
import numpy as np
import os
import csv
import requests

def download_pdb(pdb_id):
    """
    Downloads the PDB file from the RCSB PDB website if not available locally.
    
    Args:
        pdb_id (str): The PDB ID to download.
    
    Returns:
        str: The filename of the downloaded PDB file.
    """
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    pdb_file = f"pdbs/{pdb_id}.pdb"
    
    if not os.path.exists(pdb_file):
        print(f"Downloading {pdb_id} from RCSB PDB...")
        response = requests.get(pdb_url)
        if response.status_code == 200:
            with open(pdb_file, 'w') as f:
                f.write(response.text)
            print(f"Downloaded {pdb_id}.pdb")
        else:
            print(f"Error: Could not download PDB file {pdb_id}.")
            return None
    return pdb_file

def calculate_distances(pdb_id, cutoff_distance=None):
    """
    Calculates distances between Trp-indole N and Tyr-O phenolic/Met S residues.
    
    Args:
        pdb_id (str): The PDB ID of the structure.
        cutoff_distance (float, optional): Maximum distance (in Å) to select Tyr/Met residues
                                           near Trp-indole N. If None, select all.
    Returns:
        tuple: Lists of tuples with distances to Tyr and Met (e.g., [('1FTN', 'W99', 'Y160', 6.43)]).
    """
    pdb_file = download_pdb(pdb_id)
    if pdb_file is None:
        return [], []  # Skip if download failed
    
    u = mda.Universe(pdb_file)
    
    # Select Trp-indole N
    trp_n = u.select_atoms("resname TRP and name N")
    
    tyr_o = u.select_atoms("resname TYR and name OH")
    met_s = u.select_atoms("resname MET and name SD")

    tyr_results = []
    met_results = []

    if len(trp_n) == 0:
        print(f"No Trp-indole N atoms found in {pdb_id}.")
        return tyr_results, met_results

    # Calculate distances for Trp to Tyr
    for trp_atom in trp_n:
        for tyr_atom in tyr_o:
            distance = np.linalg.norm(trp_atom.position - tyr_atom.position)
            if cutoff_distance is None or distance <= cutoff_distance:
                trp_resid = f"W{trp_atom.resid}"
                tyr_resid = f"Y{tyr_atom.resid}"
                tyr_results.append((pdb_id, trp_resid, tyr_resid, distance))

    # Calculate distances for Trp to Met
    for trp_atom in trp_n:
        for met_atom in met_s:
            distance = np.linalg.norm(trp_atom.position - met_atom.position)
            if cutoff_distance is None or distance <= cutoff_distance:
                trp_resid = f"W{trp_atom.resid}"
                met_resid = f"M{met_atom.resid}"
                met_results.append((pdb_id, trp_resid, met_resid, distance))

    return tyr_results, met_results

def process_pdb_list(pdb_list_file, tyr_csv_output_file, met_csv_output_file, cutoff_distance=None):
    """
    Reads a file containing PDB IDs and calculates distances for each PDB.
    
    Args:
        pdb_list_file (str): Path to the text file containing the PDB IDs (one per line).
        tyr_csv_output_file (str): Path to the CSV file for Tyr distances.
        met_csv_output_file (str): Path to the CSV file for Met distances.
        cutoff_distance (float, optional): Maximum distance (in Å) to select Tyr/Met residues
                                           near Trp-indole N. If None, select all.
    """
    with open(pdb_list_file, 'r') as f:
        pdb_ids = [line.strip() for line in f.readlines()]

    with open(tyr_csv_output_file, mode='w', newline='') as tyr_csvfile, \
         open(met_csv_output_file, mode='w', newline='') as met_csvfile:
        tyr_csv_writer = csv.writer(tyr_csvfile)
        met_csv_writer = csv.writer(met_csvfile)
        # Write header
        tyr_csv_writer.writerow(["PDB ID", "TRP Residue", "TYR Residue", "Distance (Å)"])
        met_csv_writer.writerow(["PDB ID", "TRP Residue", "MET Residue", "Distance (Å)"])

        for pdb_id in pdb_ids:
            print(f"Processing PDB ID: {pdb_id}")
            tyr_distances, met_distances = calculate_distances(pdb_id, cutoff_distance)
            for entry in tyr_distances:
                tyr_csv_writer.writerow(entry)
            for entry in met_distances:
                met_csv_writer.writerow(entry)

# Example usage:
families = ['rho', 'rac', 'cdc']
for fam in families:
    pdb_list_file = f'{fam}.txt'  # Text file with PDB IDs (one per line)
    tyr_csv_output_file = f'output_trp_to_tyr_{fam}.csv'  # Output CSV for Tyr distances
    met_csv_output_file = f'output_trp_to_met_{fam}.csv'  # Output CSV for Met distances
    cutoff_distance = 10.0  # Optional cutoff distance, set to None for all distances
    process_pdb_list(pdb_list_file, tyr_csv_output_file, met_csv_output_file, cutoff_distance)
