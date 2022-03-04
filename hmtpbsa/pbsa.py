import os
import sys


import shutil
import logging
import argparse
from lib.parameters import generate_input_file
from utils import obtain_id_from_index

'''
1. one protein file and many other ligands file
2. one complex file and a trajectory file.
'''

def set_amber_home(proc):
    cmd = 'which %s '%proc
    f = os.popen(cmd)
    text = f.read().strip()
    if not text:
        raise Exception("Command not found: %s "%proc)
    bindir = os.path.split(text)[0]
    amberhome = os.path.split(bindir)[0]
    return amberhome

gmx_MMPBSA='gmx_MMPBSA'
if 'AMBERHOME' not in os.environ:
    #os.environ['AMBERHOME'] = set_amber_home(gmx_MMPBSA)
    raise Exception("Not found variable AMBERHOME")

def obtain_num_of_frame(trajfile):
    """
    Get the number of frames in a trajectory file
    
    Args:
      trajfile: the trajectory file, in .xtc or .trr format.
    
    Returns:
      the number of frames in the trajectory file.
    """
    cmd = 'gmx check -f %s 2>&1 |grep Coords'%trajfile
    fr = os.popen(cmd)
    text = fr.read().strip()
    if not text:
        print(cmd)
        raise Exception("ERROR obtain %s's frame number.")
    nframe = int(text.split()[1])
    return nframe
    


class PBSA(object):
    def __init__(self, workdir='MMPBSA', mode='gb') -> None:
        self.workdir = os.path.abspath(workdir)
        self.mode = mode
        self.cwd = os.getcwd()

    def set_paras(self, complexfile, trajectoryfile, topolfile, indexfile, decompose=False):
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir)
        complexfile = os.path.abspath(complexfile)
        indexfile = os.path.abspath(indexfile)
        trajectoryfile = os.path.abspath(trajectoryfile)
        topolfile = os.path.abspath(topolfile)
        os.chdir(self.workdir)
        nframe = obtain_num_of_frame(trajectoryfile)
        mmbpsafile = generate_input_file(self.mode, decompose=decompose, outfile='mmpbsa.in', startFrame=1, endFrame=nframe, interval=1, temperature=300, igbValue=2, name='Calculate')
        # mode='gb', outfile='mmpbsa.in', startFrame=1, endFrame=1, interval=1, temperature=300, igbValue=2, name='Calculate', decompose=False
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
        #print("="*80)
        logging.info('Run the gmx_MMPBSA: %s'%self.mode)
        cmd = '{gmx_MMPBSA} -i {mmpbsa} -cs {complexfile} -ci {indexfile} -ct {trajectoryfile} -cp {topolfile} -cg {receptor} {ligand} -nogui >mmpbsa.log 2>&1 '.format(**self.paras)
        RC = os.system(cmd)
        if RC != 0:
            raise Exception('ERROR run: %s \nPlease ckeck the log file for details: mmpbsa.log'%cmd)
        shutil.copy('FINAL_RESULTS_MMPBSA.dat', self.cwd)
        outfile = "FINAL_RESULTS_MMPBSA.dat"
        self.clean(verbose=verbose)
        #print('='*80)
        print('Results: %s'%outfile)

    def clean(self, verbose=0):
        #print('='*80)
        logging.info('Clean the results.')
        if not verbose:
            cmd = '{gmx_MMPBSA} --clean >>mmpbsa.log 2>&1 '.format(gmx_MMPBSA=gmx_MMPBSA)
            RC = os.system(cmd)
            shutil.rmtree(self.workdir)
        os.chdir(self.cwd)

    def extract_result(self, energyfile='FINAL_RESULTS_MMPBSA.dat'):
        tagName = False#["GENERALIZED BORN", "POISSON BOLTZMANN"]
        tagG = False
        detal_G = {}

        with open(energyfile) as fr:
            for line in fr:
                if line.startswith('GENERALIZED BORN'):
                    tagName = 'GB'
                elif line.startswith('POISSON BOLTZMANN'):
                    tagName = 'PB'
                elif line.startswith('DELTA TOTAL'):
                    if tagName:
                        lineTemp = line.split()
                        detal_G[tagName] = (float(lineTemp[2]), float(lineTemp[3]))
                    else:
                        logging.warning("Found a DELTA G without name!")
        return detal_G

def main():
    parser = argparse.ArgumentParser(description='Free energy calcaulated by PBSA method.')
    parser.add_argument('-i', dest='INP', help='A pdb file or a tpr file for the trajectory.', required=True)
    parser.add_argument('-p', dest='TOP', help='Gromacs topol file for the system.', required=True)
    parser.add_argument('-ndx', dest='ndx', help='Gromacs index file, must contain recepror and ligand group.', required=True)
    parser.add_argument('-o', dest='OUTP', help='Output floder to save results.', default=None)
    parser.add_argument('-m', dest='MODE', help='Method to calculate: gb, pb, pb+gb. default:gb', default='gb', choices=['gb', 'pb', 'pb+gb', 'gb+pb'])
    parser.add_argument('-t', dest='TRAJ', help='A trajectory file contains many structure frames. File format: xtc, pdb, gro...', default=None)
    parser.add_argument('-dec', dest='dec', help='Decompose the energy. default:false', action='store_true', default=False)
    parser.add_argument('-D', dest='DEBUG', help='DEBUG model, keep all the files.', default=False, action='store_true')

    args = parser.parse_args()
    #exit(0)
    complexFile, topolFile, indexFile,  outdir, trajFile, debug, dec, mode = args.INP, args.TOP, args.ndx, args.OUTP, args.TRAJ, args.DEBUG, args.dec, args.MODE
    if trajFile is None:
        trajFile = complexFile
    #if len(sys.argv[1:])!=5:
    #    print('Usage: python %s <complexfile> <trajfile> <topfile> <indexfile> <mode>'%sys.argv[1])
    #    exit(0)
    #complexFile, trajFile, topolFile, indexFile, mode= sys.argv[1:]
    pbsa = PBSA(mode=mode)
    pbsa.set_paras(complexfile=complexFile, trajectoryfile=trajFile, topolfile=topolFile, indexfile=indexFile, decompose=dec)
    pbsa.run(verbose=debug)
    detal_G = pbsa.extract_result()
    print("mode    detal_G(kcal/mole)    Std. Dev.")
    for k, v in detal_G.items():
        print('%4s    %18.4f    %9.4f'%(k, v[0], v[1]))


if __name__ == "__main__":
    main()