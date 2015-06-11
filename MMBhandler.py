import os, sys
import shutil
import subprocess

LD_LIBRARY_PATH = '/home/rcc/MMB2_13'

def writeCommandsdatforMMB(commandsname,position,residueToMutate,around='2'):
    fout = open(commandsname,'w')
    fout.write('loadSequencesFromPdb\n')
    fout.write('substituteResidue A '+ str(position) +' ' + residueToMutate)
    commandsbody = '''
removeRigidBodyMomentum TRUE
mobilizer Rigid
'''+'mobilizer Default A '+position+'-'+around+' '+ position+'+'+around+'''

#constrainChainRigidSegments

firstStage 2 

lastStage 2
smallGroupInertiaMultiplier 11  
reportingInterval .1
numReportingIntervals 150 
useOpenMMAcceleration true

readAtStage 2
setDefaultMDParameters

readBlockEnd

randomizeInitialVelocities false
'''+'includeAllResiduesWithin 1.2 A '+position+'-'+around +'\nincludeAllResiduesWithin 1.2 A '+position +'\nincludeAllResiduesWithin 1.2 A '+position+'-'+around+'''


readFromStage 1000    

readBlockEnd   
    '''
    fout.write('%s\n' % commandsbody)
    fout.close()
    
def mutatePositionWithMMB(PDBname,position,residueToMutate):
    print 'export LD_LIBRARY_PATH=%s' % LD_LIBRARY_PATH
    os.system('export LD_LIBRARY_PATH=' + LD_LIBRARY_PATH )
    os.environ['LD_LIBRARY_PATH'] = LD_LIBRARY_PATH
    os.putenv("LD_LIBRARY_PATH",LD_LIBRARY_PATH)
    os.system("export LD_LIBRARY_PATH")
    cmdsname = 'tmpcommand.dat'
    writeCommandsdatforMMB(cmdsname,position,residueToMutate,around='2')
    comando = os.path.join(LD_LIBRARY_PATH,'MMB.2_13.Linux64') + ' -c ' + cmdsname
    print comando
    shutil.copy(PDBname,'last.1.pdb')
    subprocess.check_call(comando, env=os.environ,shell=True)
    shutil.copy('last.2.pdb',PDBname[:-4]+position+residueToMutate+'.pdb')
    
def mutatePositionsWithMMB(PDBname,position):
    aminolist = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    for aa in aminolist:
        mutatePositionWithMMB(PDBname,position,aa)
    
if __name__ == '__main__':
    PDBname,position = sys.argv[1], str(sys.argv[2])
    mutatePositionsWithMMB(PDBname,position)
