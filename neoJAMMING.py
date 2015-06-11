from prody import *
from collections import defaultdict
import RCCpackage.RCCutils as rccu
import RCCpackage.RCCdata as rccd
from pymin import *
import sys
import operator
import os
import shutil
import os.path

def combineChains(pdbslist=None, chainslist=None,newpdbname='combinedPDBs.pdb',newchain='A'):
    fout = open(newpdbname,'w')
    newresnum = -1
    atomnum = 0
    originalnames = dict()
    for pdb, chain in zip(pdbslist,chainslist):
        lastres = -1

        fin = open(pdb,'r')
        for l in fin:
            if not l.startswith('ATOM'): continue
            _chain = l[21]
            if _chain != chain: continue
            l = l.strip()
            ll = l
            #print l
            currentres = int(l[22:26])
            if lastres < currentres:
                newresnum += 1
            print 'newresnum', newresnum
            lastres = currentres
            aa3 = l[17:20]
            aa1 = rccd.AA_three_to_one.get(aa3,'?')
            atomnum += 1
            satomnum = '% *d' % (7, atomnum) # ('{0:07d}'.format(atomnum)).replace('0',' ')
            sresnum = '% *d' % (4, newresnum) # ('{0:04d}'.format(newresnum)).replace('0',' ')
            #ll = l[:4] + satomnum + l[11:17]+ '-' + aa1 + chain + l[20] + newchain + sresnum + l[26:]
            ll = l[:4] + satomnum + l[11:17]+ aa3 + l[20] + newchain + sresnum + l[26:]
            originalname = (aa3+str(currentres)).replace(' ','')
            originalnames[(aa3 + sresnum).replace(' ','')] = (originalname,_chain)
            fout.write('%s\n' % ll)
        fin.close()
    fout.write('END')
    fout.close()
    return originalnames

def neoJAMMING(pdbname,chain,nummodes = 10, numconfs=1000,rmsd=3.8):
    prot = parsePDB(pdbname,chain=chain)
    anm = ANM(pdbname + ' ANM analysis')
    calphas = prot.select('calpha')
    anm.buildHessian(calphas)
    anm.calcModes(nummodes)
    nmdfn = pdbname+'.nmd'
    bb_anm, bb_atoms = extendModel(anm, calphas, prot.select('protein'))
    #writeNMD(nmdfn, bb_anm, bb_atoms)
    writeNMD(nmdfn, bb_anm, bb_atoms ) #anm[:3], calphas)
    ensemble = sampleModes(bb_anm[:3], bb_atoms, n_confs=numconfs, rmsd=rmsd)
    backbone = bb_atoms.copy()
    backbone.addCoordset(ensemble)
    pdbensamblefn = nmdfn+'.pdb'
    writePDB(pdbensamblefn, backbone)
    #rccu.disassemble_NMR(pdbensamblefn)

    _basename = pdbname[:pdbname.find('.')] #individual pdbs are in this directory
    if os.path.exists(_basename):
        pass
        #shutil.rmtree(_basename)
    else:
        os.mkdir(_basename)

        for i in range(1, backbone.numCoordsets()):  # skipping 0th coordinate set
            fn = os.path.join(_basename,  str(i) + os.path.basename(pdbname))
            writePDB(fn, backbone, csets=i)

    residuescores = defaultdict(float)
    individualResScores = []
    for f in rccu.iter_directory_files(_basename):
        tmpresiduescores = defaultdict(float)
        print f , '**********',
        it = neoJAMMINGresidueScoresIter(f, chain)
        print '*******',
        for i in it:
            residuescores[i[0]] += i[1]
            tmpresiduescores[i[0]] = i[1]
        individualResScores.append(sorted(tmpresiduescores.iteritems(), key=operator.itemgetter(1),reverse=True))
        print '*****'

    return dict(all=sorted(residuescores.iteritems(), key=operator.itemgetter(1),reverse=True),individual=individualResScores)

def singlePDB():
    pdbfn = sys.argv[1]
    chain = sys.argv[2]
    sortedresidueScores = neoJAMMING(pdbfn,chain)['all']
    for rs in sortedresidueScores:
        print rs[0], '%.3f' % rs[1]
    print '[done]'

def multiplePDB():
    pdbschainslist = sys.argv[1:]
    pdbs, chains = [], []
    for i in xrange(0,len(pdbschainslist)-1,2):
        pdbs.append(pdbschainslist[i])
        chains.append(pdbschainslist[i+1])

    print pdbs
    print chains
    newpdbname='combinedPDBs.pdb'
    newchain = 'A'
    originalnames = combineChains(pdbslist=pdbs, chainslist=chains,newpdbname=newpdbname,newchain=newchain)

    sortedresidueScores = neoJAMMING(newpdbname,newchain)['all']
    for rs in sortedresidueScores:
        resname, chainname = originalnames[rs[0]]
        print rs[0].rjust(7) + ' '+ ('%.3f' % rs[1]).ljust(10) + resname.rjust(7) + '  ' + chainname.ljust(3)
    print 'done'


if __name__== '__main__':
    singlePDB()
    #multiplePDB()
