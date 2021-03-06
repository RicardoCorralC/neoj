#python permutationDistanceKDE.py 1POT.pdb A 1A99.pdb  A
import RCCpackage.RCCutils as rcu
import neoJAMMING as neoj
import sys
from scipy.spatial.distance import minkowski
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from sklearn import manifold
import numpy as np
from sklearn.manifold import Isomap
from sklearn.manifold import spectral_embedding
from sklearn.manifold import LocallyLinearEmbedding
from sklearn.manifold import MDS
from sklearn.random_projection import SparseRandomProjection
from scipy.stats import spearmanr, pearsonr, ranksums, mannwhitneyu

from scipy import stats

import Hausdorff as hausdorff

def get_permutation(jamminglist, aligned_residues):
    #print 'jamminglist', jamminglist
    #print 'aligned_residues', aligned_residues
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
    
def get_permutations(jamminglists, aligned_residues):
    #print 'jamminglist', jamminglist
    #print 'aligned_residues', aligned_residues
    std_perums = []
    for jamminglist in jamminglists:
        permu = []
        std_permu = [0]*len(aligned_residues)
        j = 0
        for res in jamminglist:
            if res[0] not in aligned_residues: continue
            i_ = aligned_residues.index(res[0])
            permu.append(i_)
            std_permu[i_] = j
            j += 1
        std_perums.append(std_permu)
    return std_perums
    
    
def scoreFromPvalues(d1,d2):
    '''
    Here, a nonparametric test is made over each d1[i], d2[i] pair of components to quantify the overall similarity of distributions.
    '''
    pvals = 0
    n = len(d1[0]) #This should be the core size
    print 'n, the core size: ',   n
    for i in xrange(n):
        print d1[:,i],d2[:,i]
        pvals += spearmanr(d1[:,i],d2[:,i])[0] #(Pearson's correlation coefficient,2-tailed p-value)
    pvals /= n
    return pvals
    
def scoreFromKDE(d1,d2):
    '''
    Here, a nonparametric test is made over each d1[i], d2[i] pair of components to quantify the overall similarity of distributions.
    '''
    kernel = stats.gaussian_kde(d1.T, bw_method='silverman')
    print 'kpdf dim, n', kernel.d, kernel.n
    loglike = 0
    n = len(d1[0]) #This should be the core size
    print 'n, the core size: ',   n
    for p in d2:
        print p
        print 'Kernel(p)', kernel(p) 
        loglike += np.log(kernel(p)) 
    return loglike,kernel

if __name__ == '__main__':
    pnorm = 1
    pdb1, chain1 = sys.argv[1], sys.argv[2]
    pdb2, chain2 = sys.argv[3], sys.argv[4]
    nj1 = neoj.neoJAMMING(pdb1,chain1,nummodes=5,numconfs=1000,rmsd=3.8)
    nj2 = neoj.neoJAMMING(pdb2,chain2,nummodes=5,numconfs=1000,rmsd=3.8)
    neojamming1 = nj1['all']
    neojamming2 = nj2['all']

    
    offset1 = rcu.pdb_residue_offset(pdb1,chain1)
    offset2 = rcu.pdb_residue_offset(pdb2,chain1)
    print 'offset1: '  , offset1 
    print 'offset2: '  , offset2 
    
    #HERE structures must have only atoms of selected chain
    TM_align = rcu.TM_aligned_residues(pdb1,pdb2,offset1, offset2)
    
    
    individualjammings1 = np.asarray(get_permutations(nj1['individual'],TM_align['alignedList1']))
    individualjammings2 = np.asarray(get_permutations(nj2['individual'],TM_align['alignedList2']))
    

        
    clf = LocallyLinearEmbedding(n_components=2)#,method='modified')#Isomap(n_components=2)
    clf.fit(individualjammings1)
    ij1 = clf.transform(individualjammings1)
    ij2 = clf.transform(individualjammings2)
    PValsScore,KDE = scoreFromKDE(ij1,ij2)
    print 'PValsScore: ', PValsScore
    H = hausdorff.H(individualjammings1,individualjammings2)
    print 'Hausdorff distance', H
    
    print ij1
    f, (ax1, ax2,ax3) = pl.subplots(1,3, sharex=True, sharey=True)
    pl.ioff()
    #pl.axis('off')
    pl.tight_layout(pad=0)
    pl.title('log-likelihood: %.4f\tHausdorff: %.4f'%(PValsScore,H))
    #pl.subplot(1,2,1)

    xmin = ij1[:,0].min()
    xmax = ij1[:,0].max()
    ymin = ij1[:,1].min()
    ymax = ij1[:,1].max()
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    Z = np.reshape(KDE(positions).T, X.shape)
    ax1.imshow(np.rot90(Z), cmap=pl.cm.Spectral_r,aspect=1,extent=[xmin, xmax, ymin, ymax])
    ax1.scatter(ij1[:,0],ij1[:,1],marker='o',s=2,facecolor='0.4',edgecolor='r')

    #pl.subplot(1,2,2)
    ax2.scatter(ij2[:,0],ij2[:,1],marker='o',s=2,facecolor='0.4',edgecolor='r')
    ax3.scatter(ij2[:,0],ij2[:,1],marker='o',s=12,facecolor='y',edgecolor='0.05',alpha=0.6)
    ax3.scatter(ij1[:,0],ij1[:,1],marker='o',s=12,facecolor='b',edgecolor='0.05',alpha=0.5)
    ax1.axes.get_xaxis().set_visible(False)
    ax2.axes.get_xaxis().set_visible(False)
    ax3.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.set_xlim([xmin, xmax])
    ax1.set_ylim([ymin, ymax])
    ax2.set_xlim([xmin, xmax])
    ax2.set_ylim([ymin, ymax])
    ax3.set_xlim([xmin, xmax])
    ax3.set_ylim([ymin, ymax])
    ax1.axis('off')
    #ax2.axis('off')
    #ax3.axis('off')
    pl.savefig(pdb1+'vs'+pdb2+'5mLLE.png', bbox_inches='tight',pad_inches=0, dpi=250)
    pl.show()
    pl.clf()
    
    print TM_align['seqA']
    print TM_align['matchs']
    print TM_align['alignedList1']
    print 'Aligned positions1:\n', ','.join(i[3:] for i in TM_align['alignedList1'])
    print TM_align['alignedList2']
    print 'Aligned positions2:\n', ','.join(i[3:] for i in TM_align['alignedList2'])

    
    permu1,stdpermu1 = get_permutation(neojamming1, TM_align['alignedList1'])
    permu2,stdpermu2 = get_permutation(neojamming2, TM_align['alignedList2'])
    pl.title('Conformer ensemble distance: %.2f'%minkowski(stdpermu1,stdpermu2,pnorm))
    pl.scatter(stdpermu1,stdpermu2,marker='o',s=55,facecolor='0.6',edgecolor='b')
    pl.xlim(0,len(stdpermu1))
    pl.ylim(0,len(stdpermu1))
    pl.savefig(pdb1+'vs'+pdb2+'.png', bbox_inches='tight',dpi=270)
    pl.show()
    
    
    print 'permu:\n',','.join([str(i) for i in permu1])
    print 'std permu:\n',','.join([str(i) for i in stdpermu1])
    print 'permu:\n',','.join([str(i) for i in permu2])
    print 'std permu:\n',','.join([str(i) for i in stdpermu2])
    
    print 'Conformation ensambles at distance:\t', minkowski(stdpermu1,stdpermu2,pnorm)
    fpermus = open(pdb1+'vs'+pdb2+'.txt','w')
    fpermus.write('%s\n%s\n' % (pdb1,','.join(map(str,stdpermu1))) )
    fpermus.write('%s\n%s\n' % (pdb2,','.join(map(str,stdpermu2))) )
    print 'Spearman rank-order correlation coefficient and p-value', spearmanr(permu1,permu2)
    
