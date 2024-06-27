from Bio.PDB import PDBParser
import numpy as np
import os

def get_sulfur_coordinates(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('peptide', pdb_file)
    
    sulfur_coords = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == 'CYS':
                    for atom in residue:
                        if atom.get_id() == 'SG':
                            sulfur_coords.append(atom.get_coord())
    
    if len(sulfur_coords) != 2:
        raise ValueError("The PDB file does not contain exactly two Cysteine residues.")
    
    return sulfur_coords

def calculate_distance(coord1, coord2):
    return np.linalg.norm(coord1 - coord2)

def main(pdb_files_directory):
    for pdb_file in os.listdir(pdb_files_directory):
        if pdb_file.endswith(".pdb"):
            file_path = os.path.join(pdb_files_directory, pdb_file)
            sulfur_coords = get_sulfur_coordinates(file_path)
            distance = calculate_distance(sulfur_coords[0], sulfur_coords[1])
            print(f"Distance between thiol groups in {pdb_file}: {distance:.2f} Å")

# Example usage
pdb_files_directory = "../data/Peptides/spr_binders/CNLWEFECGASG/"
main(pdb_files_directory)