from Bio.PDB import PDBParser
import numpy as np
import os
import argparse
import sys

parser = argparse.ArgumentParser(description="Calculate distance between sufur atoms and output bond angles",
                                 usage="python thiol_distance.py --pdb_dir directory_path")
parser.add_argument("--pdb_dir", help="folder containing PDB files")
parser.add_argument("--output", help="Output file", default="output.tsv")
args = parser.parse_args()

def get_cysteine_atoms(pdb_file: str) -> list:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('peptide', pdb_file)
    # we only have one model and one chain
    chain = structure[0]["A"]
   
    cys_atoms = []
        
    chain.atom_to_internal_coordinates(verbose=True)
    for residue in chain:
        if residue.get_resname() == 'CYS':            
            ca_atom = residue['CA']
            cb_atom = residue['CB']
            sg_atom = residue['SG']
            cys_atoms.append((ca_atom.get_coord(), cb_atom.get_coord(), sg_atom.get_coord()))

    if len(cys_atoms) != 2:
        raise ValueError("The PDB file does not contain exactly two Cysteine residues.")
    
    return cys_atoms

def calculate_distance(coord1: np.ndarray, coord2: np.ndarray) -> float:
    # could also use
    # from sklearn.metrics.pairwise import paired_euclidean_distances
    # paired_euclidean_distances(coord1,coord2)[0]
    return np.linalg.norm(coord1 - coord2)

# next three functions are from https://pycrawfordprogproj.readthedocs.io/en/latest/Project_01/Project_01.html
def bond_unit_vector(i: np.ndarray, j: np.ndarray) -> np.ndarray:
    # Input: `i`, `j` index of molecule's atom
    # Output: Unit vector of bond from atom `i` to atom `j`
    vec = j - i
    return vec / np.linalg.norm(vec)

def bond_angle(i: np.ndarray, j: np.ndarray, k: np.ndarray) -> float:
    # Input: `i`, `j`, `k` index of molecule's atom; where `j` is the central atom
    # Output: Bond angle for atoms `i`-`j`-`k`
    e_ji = bond_unit_vector(j, i)
    e_jk = bond_unit_vector(j, k)
    return np.arccos(e_ji.dot(e_jk)) * 180 / np.pi

def dihedral_angle(i: np.ndarray, j: np.ndarray, k: np.ndarray, l: np.ndarray) -> float:
    # Input: `i`, `j`, `k`, `l` index of molecule's atom; where `k` is the central atom, and angle is i - j-k-l
    # Output: Dihedral angle for atoms `i`-`j`-`k`-`l`
    # This alternative code uses cross products of the bond unit vectors to find the normal vectors to the planes formed by the atoms, then uses the dot product of these normal vectors.
    # It uses cross products of the bond unit vectors
    # The alternative method explicitly computes bond angles and uses these to adjust the dihedral angle calculation.
    res = np.cross(bond_unit_vector(j, i), bond_unit_vector(j, k)).dot(np.cross(bond_unit_vector(k, j), bond_unit_vector(k, l)))
    res /= np.sin(bond_angle(i, j, k) / 180 * np.pi) * np.sin(bond_angle(j, k, l) / 180 * np.pi)
    assert(np.abs(res) < 1 + 1e-7)
    res = np.sign(res) if np.abs(res) > 1 else res
    return np.arccos(res) * 180 / np.pi


    
def process_pdbs(pdb_files_directory: str) -> list:

    outstrings = []
    for pdb_file in os.listdir(pdb_files_directory):
        if pdb_file.endswith(".pdb"):
            file_path = os.path.join(pdb_files_directory, pdb_file)
            try:
                cys_atoms = get_cysteine_atoms(file_path)
                # Coordinates of CA and SG atoms
                ca1, cb1, sg1 = cys_atoms[0]
                ca2, cb2, sg2 = cys_atoms[1]
                
                # Calculate distance between SG atoms
                distance = calculate_distance(sg1, sg2)
                outstrings.append("Distance between thiol groups in {}: {:.2f} Å".format(pdb_file,distance))
                
                # Calculate bond angle CA-SG-CA
                angle_ca1_sg1_ca2 = bond_angle(ca1, sg1, ca2)
                angle_ca1_sg2_ca2 = bond_angle(ca1, sg2, ca2)
                outstrings.append("Angle between CA1-SG1-CA2 in {}: {:.2f} degrees".format(pdb_file,angle_ca1_sg1_ca2))
                outstrings.append("Angle between CA1-SG2-CA2 in {}: {:.2f} degrees".format(pdb_file,angle_ca1_sg2_ca2))
                # Calculate bond angle CB-SG-CB
                angle_cb1_sg1_cb2 = bond_angle(cb1, sg1, cb2)
                angle_cb1_sg2_cb2 = bond_angle(cb1, sg2, cb2)
                outstrings.append("Angle between CB1-SG1-CB2 in {}: {:.2f} degrees".format(pdb_file,angle_cb1_sg1_cb2))
                outstrings.append("Angle between CB1-SG2-CB2 in {}: {:.2f} degrees".format(pdb_file,angle_cb1_sg2_cb2))
                # Calculate bond angle CA-CB-SG
                angle_ca1_cb1_sg1 = bond_angle(ca1, cb1, sg1)
                angle_ca2_cb2_sg2 = bond_angle(ca2, cb2, sg2)                
                outstrings.append("Angle between CA1-CB1-SG1 in {}: {:.2f} degrees".format(pdb_file,angle_ca1_cb1_sg1))
                outstrings.append("Angle between CA2-CB2-SG2 in {}: {:.2f} degrees".format(pdb_file,angle_ca2_cb2_sg2))

                # Calculate dihedral angle CB-SG-SG-CB
                dihedral_angle_ca_sg_sg_ca = dihedral_angle(ca1, sg1, sg2, ca2)
                dihedral_angle_cb_sg_sg_cb = dihedral_angle(cb1, sg1, sg2, cb2)                
                outstrings.append("Dihedral angle CA-SG-SG-CA in {}: {:.2f} degrees".format(pdb_file,dihedral_angle_ca_sg_sg_ca))
                outstrings.append("Dihedral angle CB-SG-SG-CB in {}: {:.2f} degrees\n".format(pdb_file,dihedral_angle_cb_sg_sg_cb))
            
            except ValueError as e:
                print(e)

    return outstrings

def main():

    pdb_files_directory = args.pdb_dir # "../data/Peptides/spr_binders/ABDCGG/"
    bond_strings = process_pdbs(pdb_files_directory)

    try:
        with open(args.output, "w") as f:
            f.write('\n'.join(bond_strings))
    except IOError as e:
        print(f"Error writing to file {filename}: {e}")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
