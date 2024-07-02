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
            sg_atom = residue['SG']
            ca_atom = residue['CA']
            cys_atoms.append((ca_atom.get_coord(), sg_atom.get_coord()))
            #if atom.get_id() == 'SG':
            #    sulfur_coords.append(atom.get_coord())

    if len(cys_atoms) != 2:
        raise ValueError("The PDB file does not contain exactly two Cysteine residues.")
    
    return cys_atoms

def calculate_distance(coord1: np.ndarray, coord2: np.ndarray) -> float:
    # could also use
    # from sklearn.metrics.pairwise import paired_euclidean_distances
    # paired_euclidean_distances(coord1,coord2)[0]
    return np.linalg.norm(coord1 - coord2)

def calculate_angle(coord1: np.ndarray, coord2: np.ndarray, coord3: np.ndarray) -> float:
    vec1 = coord1 - coord2
    vec2 = coord3 - coord2
    cosine_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)

def calculate_dihedral(coord1: np.ndarray, coord2: np.ndarray, coord3: np.ndarray, coord4: np.ndarray) -> float:
    # This function uses vector projections and cross products to compute the dihedral angle.
    # It projects vectors onto a plane orthogonal to one of the bond vectors (b2), providing a clear geometric interpretation of the dihedral angle, which tends to be more widely used and straightforward for structural biology contexts

    # Vectors between points
    b1 = coord2 - coord1
    b2 = coord3 - coord2
    b3 = coord4 - coord3

    # Normalize b2 so that it does not influence magnitude of vector rejections
    b2 /= np.linalg.norm(b2)

    # Orthogonal vectors to b2
    v = b1 - np.dot(b1, b2) * b2
    w = b3 - np.dot(b3, b2) * b2

    # Normal vectors
    x = np.dot(v, w)
    y = np.dot(np.cross(b2, v), w)
    
    angle = np.degrees(np.arctan2(y, x))
    return angle

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
                ca1, sg1 = cys_atoms[0]
                ca2, sg2 = cys_atoms[1]
                
                # Calculate distance between SG atoms
                distance = calculate_distance(sg1, sg2)
                
                # Calculate bond angle CA-SG-CA
                angle_ca_sg_ca = calculate_angle(ca1, sg1, ca2)
                angle_ca_sg_ca2 = bond_angle(ca1, sg1, ca2)
                
                # Calculate dihedral angle CA-SG-SG-CA
                dihedral_angle_ca_sg_sg_ca = calculate_dihedral(ca1, sg1, sg2, ca2)
                dihedral_angle_ca_sg_sg_ca2 = dihedral_angle(ca1, sg1, sg2, ca2)                                

                outstrings.append("Distance between thiol groups in {}: {:.2f} Å".format(pdb_file,distance))
                outstrings.append("Angle between CA-SG-CA in {}: {:.2f} degrees".format(pdb_file,angle_ca_sg_ca))
                outstrings.append("Dihedral angle CA-SG-SG-CA in {}: {:.2f} degrees - vector projection onto a plane orthogonal a bond vector".format(pdb_file,dihedral_angle_ca_sg_sg_ca))
                #outstrings.append("Angle between CA-SG-CA in {}: {:.2f} degrees".format(pdb_file,angle_ca_sg_ca2))
                outstrings.append("Dihedral angle CA-SG-SG-CA in {}: {:.2f} degrees - adjusted by explicit bond angle computation\n".format(pdb_file,dihedral_angle_ca_sg_sg_ca2))
            
            except ValueError as e:
                print(e)

    return outstrings

def main():

    pdb_files_directory = args.pdb_dir # "../data/Peptides/spr_binders/CNLWEFECGASG/"
    bond_strings = process_pdbs(pdb_files_directory)

    try:
        with open(args.output, "w") as f:
            f.write('\n'.join(bond_strings))
    except IOError as e:
        print(f"Error writing to file {filename}: {e}")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())