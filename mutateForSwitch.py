#This program takes two structures A and B. 
#A target structure A is used to select the mutant structure of B with the most similar conformer population

import RCCpackage.RCCutils as rcu
import neoJAMMING as neoj
import sys
from scipy.spatial.distance import euclidean, minkowski
from scipy.stats import spearmanr
import matplotlib.pylab as pl
import random
import os



NUM_MUTANTS = 100
MUTA_PROBA = 0.5


def generateMutantSeq(S,aligned,p):
	'''
	A Bernoulli trial with probability p is carried out in each position of S not aligned.
	'''
	aminolist = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	newseq = []
	for i in xrange(len(S)):
		if S[i] == '-': continue
		newchar = S[i]
		if aligned[i] != ':':
			if random.random()<=p:
				newchar = random.choice(aminolist)
		newseq.append(newchar)
	nS = ''.join(newseq)
	print 'Original Sequence:\n', S
	print 'Mutated Sequence:\n', nS
	return nS


def generateMutantSeqList(S,aligned,p,N):
	mutatedSeqs = []
	for i in xrange(N):
		mutatedSeqs.append(generateMutantSeq(S,aligned,p))
	return mutatedSeqs


def makeRunFileForFoldX(pdbfile,S,aligned,p,N):
	'''
	Dirty function to make necessary files for FoldX
	'''
	os.system('echo \''+pdbfile+'\' > list.txt')
	mutatedSeqs = generateMutantSeqList(S,aligned,p,N)
	#os.system('rm -rf mutant_file.txt')
	#os.system('touch mutant_file.txt')
	originalSeq = S.replace('-','')
	os.system('echo \''+originalSeq+'\' > mutant_file.txt')
	for ms in mutatedSeqs:
		os.system('echo \''+ms+'\' >> mutant_file.txt')
	os.system('FoldX -runfile run.txt') #run.txt must exist already and FoldX must be in path

def get_permutation(jamminglist, aligned_residues):
	permu = []
	std_permu = [0]*len(aligned_residues)
	j = 0
	for res in jamminglist:
		if res[0] not in aligned_residues: continue
		i_ = aligned_residues.index(res[0])
		permu.append(i_)
		std_permu[i_] = j
		j += 1

	return permu,std_permu

if __name__ == '__main__':
	pdb1, chain1 = sys.argv[1], sys.argv[2]
	pdb2, chain2 = sys.argv[3], sys.argv[4]
	
	MUTA_PROBA = float(sys.argv[5])
	
	offset1 = rcu.pdb_residue_offset(pdb1,chain1)
	offset2 = rcu.pdb_residue_offset(pdb2,chain1)
	
	#HERE structures must have only atoms of selected chain
	TM_align = rcu.TM_aligned_residues(pdb1,pdb2,offset2, offset1)
	
	print TM_align['seqA']
	print TM_align['matchs']
	print TM_align['alignedList1']
	print 'Aligned positions1:\n', ','.join(i[3:] for i in TM_align['alignedList1'])
	print TM_align['alignedList2']
	print 'Aligned positions2:\n', ','.join(i[3:] for i in TM_align['alignedList2'])
	
	makeRunFileForFoldX(pdb2,TM_align['seqA'],TM_align['matchs'],MUTA_PROBA,NUM_MUTANTS)

	
	neojamming1 = neoj.neoJAMMING(pdb1,chain1)
	neojamming2 = neoj.neoJAMMING(pdb2,chain2)
	

	print 'offset1: '  , offset1 
	print 'offset2: '  , offset2 
	

	

	
	permu1,stdpermu1 = get_permutation(neojamming1, TM_align['alignedList1'])
	permu2,stdpermu2 = get_permutation(neojamming2, TM_align['alignedList2'])
	pl.title('Conformer ensemble distance: %.2f'%minkowski(stdpermu1,stdpermu2,1))
	pl.scatter(stdpermu1,stdpermu2,marker='o',s=55,facecolor='0.6',edgecolor='b')
	pl.xlim(0,len(stdpermu1))
	pl.ylim(0,len(stdpermu1))
	pl.show()
	
	
	print 'permu:\n',','.join([str(i) for i in permu1])
	print 'std permu:\n',','.join([str(i) for i in stdpermu1])
	print 'permu:\n',','.join([str(i) for i in permu2])
	print 'std permu:\n',','.join([str(i) for i in stdpermu2])
	
	print 'Conformation ensambles at distance:\t', minkowski(stdpermu1,stdpermu2,1)
	print 'Spearman rank-order correlation coefficient and p-value', spearmanr(permu1,permu2)
