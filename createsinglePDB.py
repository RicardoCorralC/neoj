import os, sys
import RCCpackage.RCCutils as rccu

finalf = open('PDBsote.pdb','w')
modelnum = 1
for f in rccu.iter_directory_files(sys.argv[1]):
    #finalf.write('MODEL\t%d\n'%modelnum)
    modelnum+=1
    ftmp = open(f,'r')
    for l in ftmp:
        if l.startswith('REMARK'): continue
     #   if l.startswith('END'): continue
        #if l.startswith('TER'): continue
        finalf.write('%s'%l)
    #finalf.write('ENDMDL\n')
    
print 'yhea'

