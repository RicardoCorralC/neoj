# File: protnet.py
# Coded by Mauricio Cruz Loya
# Version: 0.3a
# Last update: October, 2013
# Group: IRonBios (Gabriel del Rio's lab) at IFC, UNAM

from networkx import Graph
from scipy.spatial import cKDTree
from pdb import *

class ProteinStructureNetwork(Graph):
    '''
    Abstract class representing a protein structure network.
    '''
    def __init__(self):
        super(ProteinStructureNetwork, self).__init__()

class ResidueInteractionNetwork(ProteinStructureNetwork):
    '''
    Abstract class representing a residue interaction network. This is
    essentially a NetworkX Graph that has Biopython Residues as nodes.
    '''
    def __init__(self):
        super(ResidueInteractionNetwork, self).__init__()

    def build_from_residues(self, residues):
        '''
        residues: List of Biopython Residue objects to build the network from.
        '''
        raise NotImplementedError

    def _get_res_from_pdb(self, pdb, chains):
        '''
        Gets a list of residus from a Brookhaven Protein Data Bank file.

        pdb: String representing the pdb code of the file to use.
        chains: Either a string representing the protein chain to use or a tuple
        or list of strings representing multiple chains to use.
        '''
        if chains is str:
            residues = get_aa_residues(pdb, chains)
        else:
            # Assume 'chains' is an iterable list of chains.
            residues = []
            for chain in chains:
                residues.extend(get_aa_residues(pdb, chain))
        # Returns for use in subclasses. This is not meant to be used directly.
        return residues
            
class UnweightedRIN(ResidueInteractionNetwork):
    '''
    This class is meant for a residue interaction network with unweighted edges.
    All residues within a distance cutoff are connected by an edge. Two possible
    criteria are used for the distance: a) the distance between alpha carbons or
    b) the distance between any pair of atoms.
    '''
    def __init__(self):
        super(UnweightedRIN, self).__init__()
    
    def build_from_residues(self, residues, distance_cutoff=5.0, 
                            criterion='any', verbose=False):
        '''
        residues: List of Biopython residue objects.

        distance_cutoff: Distance in angstroms to use as cutoff for edge 
        inclusion between residues.

        criterion: Criterion for considering two amino acids to be within the
        distance cutoff. 
        - 'any' for edge inclusion when any atom pair of the amino acids is 
        within the cutoff (recommended). 
        - 'alpha' for only considering the distance between alpha carbons.
        '''
        self.add_nodes_from(residues)
        self.distance_cutoff = distance_cutoff

        atom_coords = []
        res_map = {}

        if criterion != 'any' and criterion != 'alpha':
            raise ValueError('Unknown criterion.')

        if criterion == 'any':
            # Build map from atom coordinates to residue and put all atom
            # coordinates into a list.
            for residue in residues:
                # Add all non-hydrogen atom coords to list and map to residue.
                for atom in residue:
                    # Skip hydrogen atoms.
                    if atom.get_name()[0] == 'H':
                        continue
                    coords = tuple(atom.get_coord())
                    atom_coords.append(coords)
                    res_map[coords] = residue
   
        else: # (This means criterion == 'alpha')
            for residue in residues:
                # Place only alpha carbons in coordinate list and residue map.
                coords = tuple(residue['CA'].get_coord())
                atom_coords.append(coords)
                res_map[coords] = residue

        tree = cKDTree(atom_coords)
        neighbors = tree.query_pairs(self.distance_cutoff)
            
        for i, j in neighbors:
            self.add_edge(res_map[atom_coords[i]], res_map[atom_coords[j]])

    def build_from_pdb(self, pdb, chains, distance_cutoff=5.0, criterion='any',
                       verbose=False):
        residues = self._get_res_from_pdb(pdb, chains)
        self.build_from_residues(residues, distance_cutoff, criterion, verbose)

    def _within_cutoff(self, res1, res2, sqr_distance_cutoff):
        '''
        res1: BioPython Residue object.

        res2: BioPython Residue object.

        sqr_distance_cutoff: Squared distance in angstroms**2 to use as cutoff
        for edge inclusion between residues.

        Calculates whether two residues have at least an atom pair within
        the specified squared cutoff distance.

        returns: True if residues are within cutoff distance, False if they
        are not.
        '''
        for atom1 in res1:
            if atom1.get_name()[0] == 'H':
                continue
            atom1coords = atom1.get_coord()
            for atom2 in res2:
                if atom2.get_name()[0] == 'H':
                    continue
                atom2coords = atom2.get_coord()
                if (get_sqr_distance(atom1coords, atom2coords) <= 
                                     sqr_distance_cutoff):
                    return True
        return False

class WeightedRIN(ResidueInteractionNetwork):
    '''
    Abstract class.
    '''
    def __init__(self):
        super(WeightedRIN, self).__init__()

class FunctionWeightedRIN(WeightedRIN):
    '''
    Creates a protein structure network with weighted edges. The weights
    correspond to a user-specifiable function of the residue pair.
    '''
    def __init__(self):
        super(FunctionWeightedRIN, self).__init__()

    def build_from_residues(self, residues, f):
        '''
        residues: List of biopython residue objects to use as nodes
        in the protein structure network.

        f: Function that takes two (biopython) residues as arguments. 
        Either returns a number that corresponds to the weight of the 
        edge between them, or returns None for no edge between them.
        '''
        self.add_nodes_from(residues)
        self.f = f
        
        for i in range(len(residues)):
            residue1 = residues[i]
            for j in range(i + 1, len(residues)):
                residue2 = residues[j]
                edge_weight = f(residue1, residue2)
                if edge_weight != None:
                    self.add_edge(residue1, residue2, weight=edge_weight)

    def build_from_pdb(self, pdb, chains, f):
        super(FunctionWeightedRIN, self).build_from_pdb(pdb, chains)
        build_from_residues(residues, f)
