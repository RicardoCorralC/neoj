from protnet import UnweightedRIN
from pdb import *
import networkx as nx

PDB = '1hiv.pdb'
#CHAIN = ['A', 'B']
CHAIN = 'A'

residues = get_aa_residues(PDB, CHAIN)

network = UnweightedRIN()
network.build_from_residues(residues, distance_cutoff=5.0)

b = nx.betweenness_centrality(network, normalized=False)

sorted_residues = sorted(b, key=b.get, reverse=True)

print 'Residue', 'Betweenness centrality'
for residue in sorted_residues:
    print get_aa_name(residue), b[residue]
