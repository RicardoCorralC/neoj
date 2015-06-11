# Coded by Mauricio Cruz Loya
# Version: 0.2a
# Last update: September, 2013
# Group: IRonBios (Gabriel del Rio's lab) at IFC, UNAM

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Residue
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings

warnings.filterwarnings('ignore', category=PDBConstructionWarning)

def get_aa_name(res):
    '''
    res: Biopython PDB Residue object representing an amino acid.

    returns: Name of residue in three letter code + residue number format (e.g. LYS23)
    '''
    return res.get_resname() + str(res.get_id()[1])

def get_aa_residues(pdb, chain):
    '''
    pdb: Protein Data Bank file.
    chain: Chain of the PDB file.

    Get the amino acids from a protein.

    returns: List of Biopython PDB Residue objects representing the amino acids
    of the specified protein.
    '''
    parser = PDBParser()
    structure = parser.get_structure('prot', pdb )
    model = structure[0]
    if chain is str:
        chain = model[chain]
        residue_list = list(chain.get_residues())
    
    else:
        residue_list = []
        for ch in chain:
            residue_list.extend(model[ch].get_residues())

    to_remove_list = []

    for res in residue_list:
        # Store non-amino acid residues in PDB in another list.
        if res.get_id()[0] != ' ':
            to_remove_list.append(res)
    
    # Remove non-amino acid residues from original list.
    for res in to_remove_list:
        residue_list.remove(res)

    return residue_list

def get_distance(coords1, coords2):
    '''
    coords1: A set of (x,y,z) coordinates.
    coords2: A set of (x,y,z) coordinates.

    From the (x,y,z) coordinates of two objects (e.g. pair of atoms), 
    calculates the distance between them.

    returns: Distance between pair of coordinates.
    '''
    import math

    distancex = coords1[0] - coords2[0]
    distancey = coords1[1] - coords2[1]
    distancez = coords1[2] - coords2[2]

    return math.sqrt((distancex*distancex) + (distancey*distancey) 
            + (distancez*distancez))

def get_sqr_distance(coords1, coords2):
    '''
    coords1: A set of (x,y,z) coordinates.
    coords2: A set of (x,y,z) coordinates.

    From the (x,y,z) coordinates of two objects (e.g. pair of atoms), 
    calculates the distance between them.

    returns: Squared distance between pair of coordinates.
    '''
    distancex = coords1[0] - coords2[0]
    distancey = coords1[1] - coords2[1]
    distancez = coords1[2] - coords2[2]

    return (distancex*distancex) + (distancey*distancey) + (distancez*distancez)

def general_center_of_mass(pdb, chain):
    '''
    pdb: Protein Data Bank file.
    chain: Chain of the PDB file.
    
    Get the General Center of Mass from a specified protein.

    returns: Tuple of x,y,z coordinates of the GCM.
    '''

    # Molecular weight of atoms in proteins.
    mw = {'C':12.017, 'N':14.0067, 'O':15.9994, 'S':32.065, 'Se':78.96}
    
    # Get a list of all heavy atoms in the protein.
    residues = get_aa_residues(pdb, chain)
    heavy_atoms = []

    for res in residues:
        for atom in res:
            # PDB files usually do not include the hydrogen atoms of a protein, so 
            # I'm assuming there's no need to filter atoms by type.
            heavy_atoms.append(atom)

    # Coordinates for GCM numerator and denominator.
    xn, yn, zn = [0.0] * 3
    d = 0 

    for atom in heavy_atoms:
        
        if atom.get_name()[:2] == 'SE':
            element = 'Se'
        else:
            element = atom.get_name()[0]
            
        mass = mw[element]
        coords = atom.get_coord()

        xn += coords[0] * mass
        yn += coords[1] * mass   
        zn += coords[2] * mass
        d += mass

    return (xn/d, yn/d, zn/d)
