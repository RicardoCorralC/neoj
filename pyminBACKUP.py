from pymin_pdb import *
from pymin_graph import *
import networkx as nx

PDB_ = '1hiv.pdb'
CHAIN_ = 'A'
DISTANCE_CUTOFF_ = 5.0


def neoJAMMINGresidueScoresIter(PDB, CHAIN, DISTANCE_CUTOFF=5.0):
	residues = get_aa_residues(PDB, CHAIN)
	network = build_unweighted_psn(residues, DISTANCE_CUTOFF)
	nodes = network.nodes()

	b = nx.betweenness_centrality(network)

	sorted_residues_b = sorted(b, key=b.get, reverse=True)

	for aminoacid in sorted_residues_b:
	    aa_name = get_aa_name(aminoacid)
	    yield  (aa_name, b[aminoacid])
    
    
if __name__ == '__main__':
	it = neoJAMMINGresidueScoresIter(PDB_, CHAIN_, DISTANCE_CUTOFF_)
	for i in it:
		print i[0], i[1]

