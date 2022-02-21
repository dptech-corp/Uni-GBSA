from cmath import log
import os
import sys

import shutil
import logging
from lib.parameters import generate_input_file
from utils import obtain_id_from_index

'''
1. one protein file and many other ligands file
2. one complex file and a trajectory file.
'''

gmx_MMPBSA='gmx_MMPBSA'
class PBSA(object):
    def __init__(self, workdir='MMPBSA', mode='gb') -> None:
        self.workdir = os.path.abspath(workdir)
        self.mode = mode
        self.cwd = os.getcwd()

    def set_paras(self, complexfile, trajectoryfile, topolfile, indexfile):
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir)
        complexfile = os.path.abspath(complexfile)
        indexfile = os.path.abspath(indexfile)
        trajectoryfile = os.path.abspath(trajectoryfile)
        topolfile = os.path.abspath(topolfile)
        os.chdir(self.workdir)
        mmbpsafile = generate_input_file(self.mode)
        receptor, ligand = obtain_id_from_index(indexfile)
        self.paras = {
            'mmpbsa': mmbpsafile,
            'gmx_MMPBSA':gmx_MMPBSA,
            'complexfile': complexfile,
            'indexfile': indexfile,
            'trajectoryfile': trajectoryfile,
            'topolfile': topolfile,
            'receptor': receptor,
            'ligand': ligand
        }

    def run(self, verbose=0):
        print("="*80)
        print('Run the gmx_MMPBSA: %s'%self.mode)
        cmd = '{gmx_MMPBSA} -i {mmpbsa} -cs {complexfile} -ci {indexfile} -ct {trajectoryfile} -cp {topolfile} -cg {receptor} {ligand} -nogui >mmpbsa.log 2>&1 '.format(**self.paras)
        RC = os.system(cmd)
        if RC != 0:
            raise Exception('ERROR run: %s \nPlease ckeck the log file for details: mmpbsa.log'%cmd)
        shutil.copy('FINAL_RESULTS_MMPBSA.dat', self.cwd)
        outfile = "FINAL_RESULTS_MMPBSA.dat"
        self.clean(verbose=verbose)
        print('='*80)
        print('Results: %s'%outfile)

    def clean(self, verbose=0):
        print('='*80)
        print('Clean the results.')
        if not verbose:
            cmd = '{gmx_MMPBSA} --clean >>mmpbsa.log 2>&1 '.format(gmx_MMPBSA=gmx_MMPBSA)
            RC = os.system(cmd)
            shutil.rmtree(self.workdir)
        os.chdir(self.cwd)




def main():
    if len(sys.argv[1:])!=5:
        print('Usage: python %s <complexfile> <trajfile> <topfile> <indexfile> <mode>'%sys.argv[1])
        exit(0)
    complexFile, trajFile, topolFile, indexFile, mode= sys.argv[1:]
    pbsa = PBSA(mode=mode)
    pbsa.set_paras(complexfile=complexFile, trajectoryfile=trajFile, topolfile=topolFile, indexfile=indexFile)
    pbsa.run(verbose=1)

if __name__ == "__main__":
    main()