import MDAnalysis as mda
import numpy as np
import os
import csv
import requests
import logging

# Set up logger
logging.basicConfig(
    filename='output.log',  # Log file name
    level=logging.INFO,     # Log level (INFO will include all levels of log messages)
    format='%(asctime)s - %(levelname)s - %(message)s',  # Log format
)

def download_pdb(pdb_id):
    """
    Downloads the PDB or CIF file from the RCSB PDB website if not available locally.
    
    Args:
        pdb_id (str): The PDB ID to download.
    
    Returns:
        str: The filename of the downloaded PDB or CIF file.
    """
    # Try downloading the PDB file first
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    pdb_file = f"pdbs/{pdb_id}.pdb"
    
    # If the PDB file is not found, attempt to download the CIF file instead
    cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    cif_file = f"pdbs/{pdb_id}.cif"

    if not os.path.exists(pdb_file) and not os.path.exists(cif_file):
        logging.info(f"Downloading {pdb_id} from RCSB PDB...")

        # Try downloading the PDB file
        response = requests.get(pdb_url)
        if response.status_code == 200:
            with open(pdb_file, 'w') as f:
                f.write(response.text)
            logging.info(f"Downloaded {pdb_id}.pdb")
            return pdb_file
        else:
            logging.error(f"Error: Could not download PDB file {pdb_id}. Trying CIF file instead.")

            # If PDB file download fails, try downloading the CIF file
            response = requests.get(cif_url)
            if response.status_code == 200:
                with open(cif_file, 'w') as f:
                    f.write(response.text)
                logging.info(f"Downloaded {pdb_id}.cif")
                return cif_file
            else:
                logging.error(f"Error: Could not download CIF file {pdb_id}.")
                return None
    elif os.path.exists(pdb_file):
        logging.info(f"PDB {pdb_id} file found locally.")
        return pdb_file
    elif os.path.exists(cif_file):
        logging.info(f"CIF {pdb_id} file found locally.")
        return cif_file

def process_pdb_first_chain(u):
    """
    For a MDAnalysis Universe object, selects only the first chain of the protein.
    If multiple chains exist, only the first chain is retained.

    Args:
        u (MDAnalysis.Universe): Universe object containing the full structure.

    Returns:
        MDAnalysis.Universe: Universe object containing either the full structure
                             (if only one chain) or just the first chain.
    """
    # Get unique chain identifiers (segids)
    chains = list(u.segments.segids)

    if len(chains) <= 1:
        # If only one chain, return the full universe
        return u
    
    # Select only the first chain
    first_chain = u.select_atoms(f"protein and segid {chains[0]}")

    if len(first_chain) == 0:
        logging.warning(f"No valid protein chains found. Returning full structure.")
        return u

    # Create a new Universe object with just the first chain
    new_u = mda.Merge(first_chain)
    new_u.add_TopologyAttr('resnames')  # Ensure residue names are included
    new_u.add_TopologyAttr('resids')    # Ensure residue IDs are included

    return new_u

def calculate_WY_distances(trp_n, tyr_o, pdb_id, cutoff_distance):
    """Calculates distances between Trp-indole N and Tyr-O phenolic residues."""
    results = []
    for trp_atom in trp_n:
        for tyr_atom in tyr_o:
            distance = np.linalg.norm(trp_atom.position - tyr_atom.position)
            if cutoff_distance is None or distance <= cutoff_distance:
                results.append((pdb_id, f"W{trp_atom.resid}", f"Y{tyr_atom.resid}", distance))
    return results

def calculate_WM_distances(trp_n, met_s, pdb_id, cutoff_distance):
    """Calculates distances between Trp-indole N and Met S residues."""
    results = []
    for trp_atom in trp_n:
        for met_atom in met_s:
            distance = np.linalg.norm(trp_atom.position - met_atom.position)
            if cutoff_distance is None or distance <= cutoff_distance:
                results.append((pdb_id, f"W{trp_atom.resid}", f"M{met_atom.resid}", distance))
    return results

def calculate_MM_distances(met_s, pdb_id, cutoff_distance):
    """Calculates non-redundant distances between Met S atoms."""
    results = []
    num_met = len(met_s)
    for i in range(num_met):
        for j in range(i + 1, num_met):  # Avoid redundant calculations
            distance = np.linalg.norm(met_s[i].position - met_s[j].position)
            if cutoff_distance is None or distance <= cutoff_distance:
                results.append((pdb_id, f"M{met_s[i].resid}", f"M{met_s[j].resid}", distance))
    return results

def calculate_distances(pdb_id, cutoff_distance=None, first_chain_only=False):
    """
    Downloads a PDB file, initializes an MDAnalysis universe, and calculates 
    distances for Trp-Tyr, Trp-Met, and Met-Met interactions.
    
    Args:
        pdb_id (str): The PDB ID of the structure.
        cutoff_distance (float, optional): Maximum distance (in Å) for consideration. If None, selects all.
        first_chain_only (bool, optional): If True, only the first chain is considered.
    
    Returns:
        dict: Dictionary containing lists of tuples for each distance type.
    """
    pdb_file = download_pdb(pdb_id)
    if pdb_file is None:
        return {"WY": [], "WM": [], "MM": []}  # Skip if download failed

    # Initialize MDAnalysis Universe
    u = mda.Universe(pdb_file)
    # Optional first chain filtering
    if first_chain_only:
        u = process_pdb_first_chain(u)
    
    # Select atoms
    trp_n = u.select_atoms("resname TRP and name N")
    tyr_o = u.select_atoms("resname TYR and name OH")
    met_s = u.select_atoms("resname MET and name SD")

    if len(trp_n) == 0:
        logging.warning(f"No Trp-indole N atoms found in {pdb_id}.")
        return {"WY": [], "WM": [], "MM": []}

    return {
        "WY": calculate_WY_distances(trp_n, tyr_o, pdb_id, cutoff_distance),
        "WM": calculate_WM_distances(trp_n, met_s, pdb_id, cutoff_distance),
        "MM": calculate_MM_distances(met_s, pdb_id, cutoff_distance),
    }

def process_pdb_list(pdb_list_file, WY_csv_output_file, WM_csv_output_file, 
                     MM_csv_output_file, cutoff_distance=None, first_chain_only=False):
    """
    Reads a file containing PDB IDs and calculates distances for each PDB.
    
    Args:
        pdb_list_file (str): Path to the text file containing the PDB IDs (one per line).
        tyr_csv_output_file (str): Path to the CSV file for Trp-Tyr distances.
        met_csv_output_file (str): Path to the CSV file for Trp-Met distances.
        MM_csv_output_file (str): Path to the CSV file for Met-Met distances.
        cutoff_distance (float, optional): Maximum distance (in Å) to select Tyr/Met residues
                                           near Trp-indole N. If None, select all.
        first_chain_only (bool, optional): If True, only the first chain is considered.
    """
    with open(pdb_list_file, 'r') as f:
        pdb_ids = [line.strip() for line in f.readlines()]

    with open(WY_csv_output_file, mode='w', newline='') as WY_csvfile, \
         open(WM_csv_output_file, mode='w', newline='') as WM_csvfile, \
         open(MM_csv_output_file, mode='w', newline='') as MM_csvfile:
        WY_csv_writer = csv.writer(WY_csvfile)
        WM_csv_writer = csv.writer(WM_csvfile)
        MM_csv_writer = csv.writer(MM_csvfile)
        # Write header
        WY_csv_writer.writerow(["PDB ID", "TRP Residue", "TYR Residue", "Distance (Å)"])
        WM_csv_writer.writerow(["PDB ID", "TRP Residue", "MET Residue", "Distance (Å)"])
        MM_csv_writer.writerow(["PDB ID", "MET Residue", "MET Residue", "Distance (Å)"])

        for pdb_id in pdb_ids:
            logging.info(f"Processing PDB ID: {pdb_id}")
            distances_dict = calculate_distances(pdb_id, cutoff_distance, first_chain_only)
            #tyr_distances, met_distances, MM_distances = calculate_MM_distances(pdb_id, cutoff_distance)
            for entry in distances_dict["WY"]:
                WY_csv_writer.writerow(entry)
            for entry in distances_dict["WM"]:
                WM_csv_writer.writerow(entry)
            for entry in distances_dict["MM"]:
                MM_csv_writer.writerow(entry)


if __name__ == "__main__":
    # Example usage:
    family = 'test' # change this to represent the set of PDB files
    pdb_list_file = 'pdb_for_met_formatted.txt'  # Text file with PDB IDs (one per line)
    
    # make a directory for the pdbs if needed
    if not os.path.exists('pdbs'):
        os.mkdir('pdbs')

    # make a directory for the family if needed
    if not os.path.exists(family):
        os.mkdir(family)
    
    # set args
    WY_csv_output_file = f'{family}/output_WY.csv'  # Output CSV for Trp-Tyr distances
    WM_csv_output_file = f'{family}/output_WM.csv'  # Output CSV for Trp-Met distances
    MM_csv_output_file = f'{family}/output_MM.csv'  # Output CSV for Met-Met distances
    cutoff_distance = 10.0      # Optional cutoff distance, set to None for all distances
    # This arg/filtering function currently does not work
    #first_chain_only = False    # Optional flag to consider only the first chain (or not)
    
    # run the main calculation
    process_pdb_list(pdb_list_file, 
                     WY_csv_output_file, 
                     WM_csv_output_file, 
                     MM_csv_output_file,
                     cutoff_distance)
                     #first_chain_only)
