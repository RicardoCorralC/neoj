#Transfer thermal factor of each atom in A to B
import os, sys
import shutil

finA, finB = open(sys.argv[1],'r'), open(sys.argv[2],'r')
fout = open(sys.argv[2]+'tmp','w')

dA = dict()

for l in finA:
    l = l.strip()
    if not l.startswith('ATOM'):
        continue
    a = l[13:26]
    b = l[61:66]
    #print a,b
    dA[a] = b
finA.close()
    
for l in finB:
    l = l.strip()
    if not l.startswith('ATOM'):
        continue
    a = l[13:26]
    if not a in dA: continue
    nl = l[:61] + dA[a] + l[66:]
    fout.write('%s\n' % nl)
    print nl
finB.close()
shutil.move(sys.argv[2]+'tmp',sys.argv[2])
fout.close()

 
print len(dA)
    
