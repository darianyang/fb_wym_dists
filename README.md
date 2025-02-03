# fb_wym_dists

Some loose scripts to calculate distances between W-Y and W-M and M-M residues for an input list of PDB IDs.

## Requirements and Installation
To use these scripts, first set up a Python environment with the following packages installed:
- MDAnalysis (`pip install MDAnalysis`) - see note below for specifics 
- requests (`pip install requests`)

### mmCIF/PDBx Note
The script will automatically download the PDB file if it doesn't already exist. If the PDB cannot be downloaded it will try to download the CIF file instead. 

In order to be able to use MDAnalysis with mmCIF/PDBx files, this feature seems to be coming soon, but you can use it now by using the code from this PR: https://github.com/MDAnalysis/mdanalysis/pull/4712

So you would have to clone MDAnalysis locally, then checkout to this PR. Then install locally via `pip install mdanalysis/package`.

Then, because the mmCIF/PDBx reader relies on gemmi, you'll need to install this as well via `pip install gemmi`.

If you prefer to avoid having to do this, you can just use MDAnalysis normally and ensure that the list of PDBs you provide is available as PDB files, or download the PDBs yourself, and if you store them in a directory named `pdbs` in the same directory that your script runs, this works too.

If there are any runtime issues check the log file (`output.log`).

## Directory
- `\cdc-rac-rho`
    - Scripts for calculating and plotting TRP-TYR and TRP-MET distances.
- `\met-met`
    - Scripts for calculating TRP-TYR, TRP-MET, and MET-MET distances.
    - These scripts supercede the ones in `\cdc-rac-rho` in a way, I would reccomend using these for all distance calcs.
    - There is also a script available here to convert a comma separated list of PDBs into the properly formatted single column (new line separated) file needed for analysis.

---
