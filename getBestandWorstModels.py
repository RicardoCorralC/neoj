#Given a target protein A and a directory with models of mutated B, select the best and worst according to permutation distance.
#time python getBestandWorstModels.py 1OYG.pdb A  2YFS.pdb A 2YFSmodels/
import RCCpackage.RCCutils as rcu
import neoJAMMING as neoj
import sys
from scipy.spatial.distance import euclidean, minkowski
from scipy.stats import spearmanr
import matplotlib.pylab as pl
import os.path

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
    NUMMODES = 50
    PNORM = 1
    pdb1, chain1 = sys.argv[1], sys.argv[2]
    pdb2, chain2 = sys.argv[3], sys.argv[4]
    directory = sys.argv[5]
    
    offset1 = rcu.pdb_residue_offset(pdb1,chain1)
    offset2 = rcu.pdb_residue_offset(pdb2,chain1)
    print 'offset1: '  , offset1 
    print 'offset2: '  , offset2 
    
    #HERE structures must have only atoms of selected chain
    TM_align = rcu.TM_aligned_residues(pdb1,pdb2,offset2, offset1)
    
    neojamming1 = neoj.neoJAMMING(pdb1,chain1)
    neojamming2 = neoj.neoJAMMING(pdb2,chain2)
    
    print TM_align['seqA']
    print TM_align['matchs']
    print TM_align['alignedList1']
    print 'Aligned positions1:\n', ','.join(i[3:] for i in TM_align['alignedList1'])
    print TM_align['alignedList2']
    print 'Aligned positions2:\n', ','.join(i[3:] for i in TM_align['alignedList2'])

    permu1,stdpermu1 = get_permutation(neojamming1, TM_align['alignedList1'])
    permu2,stdpermu2 = get_permutation(neojamming2, TM_align['alignedList2'])
    pl.title('Conformer ensemble distance: %.2f'%minkowski(stdpermu1,stdpermu2,PNORM))
    pl.scatter(stdpermu1,stdpermu2,marker='o',s=55,facecolor='0.6',edgecolor='b')
    pl.xlim(0,len(stdpermu1))
    pl.ylim(0,len(stdpermu1))
    pl.savefig(pdb1+'vs'+pdb2+'.png', bbox_inches='tight',dpi=250)

    
    
    print 'permu:\n',','.join([str(i) for i in permu1])
    print 'std permu:\n',','.join([str(i) for i in stdpermu1])
    print 'permu:\n',','.join([str(i) for i in permu2])
    print 'std permu:\n',','.join([str(i) for i in stdpermu2])
    
    print 'Conformation ensambles at distance:\t', minkowski(stdpermu1,stdpermu2,PNORM)
    print 'Spearman rank-order correlation coefficient and p-value', spearmanr(permu1,permu2)
    
    bestDistance = minkowski(stdpermu1,stdpermu2,PNORM)
    worstDistance = minkowski(stdpermu1,stdpermu2,PNORM)
    bestModel = 'original'
    worstModel = 'original'
    
    listoffiles = [i for i in rcu.iter_directory_files(directory)]
    for mutatedmodel in rcu.iter_directory_files(directory):
        if not mutatedmodel.endswith('.pdb') : continue
        if mutatedmodel.find('nmd') != -1: continue
        if mutatedmodel + '.nmd' in listoffiles: continue
    
        neojamming2 = neoj.neoJAMMING(mutatedmodel,chain2)

        permu2,stdpermu2 = get_permutation(neojamming2, TM_align['alignedList2'])
        tmpdistance = minkowski(stdpermu1,stdpermu2,PNORM)
        if tmpdistance < bestDistance:
            bestDistance = tmpdistance 
            bestModel = mutatedmodel
            print bestDistance, bestModel
        if tmpdistance > worstDistance:
            worstDistance = tmpdistance
            worstModel = mutatedmodel
        pl.clf()
        pl.title('Conformer ensemble distance: %.2f'%tmpdistance)
        pl.scatter(stdpermu1,stdpermu2,marker='o',s=55,facecolor='0.6',edgecolor='b')
        pl.xlim(0,len(stdpermu1))
        pl.ylim(0,len(stdpermu1))
        pl.savefig(pdb1+'vs'+os.path.basename(mutatedmodel)+'.png', bbox_inches='tight',dpi=250)
        
    fout = open('bestAndWorstModels'+pdb1+'vs'+pdb2+'.txt','w')
    fout.write('Best distance: %.2f\n' % bestDistance)
    fout.write('Best model: %s\n' % bestModel)
    fout.write('Worst distance: %.2f\n' % worstDistance)
    fout.write('Worst model: %s\n' % worstModel)
