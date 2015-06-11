from protnet import UnweightedRIN
from pdb import *
import networkx as nx

PDB_ = '1hiv.pdb'
CHAIN_ = 'A'
DISTANCE_CUTOFF_ = 5.0


def neoJAMMINGresidueScoresIter(PDB, CHAIN, DISTANCE_CUTOFF=5.0):
    residues = get_aa_residues(PDB, CHAIN)

    network = UnweightedRIN()
    network.build_from_residues(residues, distance_cutoff=5.0)

    b = nx.betweenness_centrality(network, normalized=False)

    sorted_residues = sorted(b, key=b.get, reverse=True)
    
    for residue in sorted_residues:
        aa_name = get_aa_name(residue)
        yield  (aa_name, b[residue] )


if __name__ == '__main__':
    it = neoJAMMINGresidueScoresIter(PDB_, CHAIN_, DISTANCE_CUTOFF_)
    for i in it:
        print i[0], i[1]

