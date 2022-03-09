import os
import shutil
import logging
import argparse
from hmtpbsa.settings import gmx_MMPBSA
from hmtpbsa.utils import obtain_id_from_index
from hmtpbsa.pbsa.utils import obtain_num_of_frame
from hmtpbsa.pbsa.parameters import generate_input_file

'''
1. one protein file and many other ligands file
2. one complex file and a trajectory file.
'''

class PBSA(object):
    def __init__(self, workdir='MMPBSA', mode='gb') -> None:
        self.workdir = os.path.abspath(workdir)
        self.mode = mode
        self.cwd = os.getcwd()

    def set_paras(self, complexfile, trajectoryfile, topolfile, indexfile, decompose=False, indi=1.0, exdi=80.0):
        """
        The function is used to set the parameters for the MMPBSA.py script
        
        Args:
          complexfile: the complex pdb file
          trajectoryfile: the trajectory file, such as md.xtc
          topolfile: the topology file of the complex
          indexfile: the index file of the system, which contains the index of the receptor and ligand.
          decompose: whether to decompose the complex into its components. Defaults to False
        """
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir)
        complexfile = os.path.abspath(complexfile)
        indexfile = os.path.abspath(indexfile)
        trajectoryfile = os.path.abspath(trajectoryfile)
        topolfile = os.path.abspath(topolfile)
        os.chdir(self.workdir)
        nframe = obtain_num_of_frame(trajectoryfile)
        mmbpsafile = generate_input_file(self.mode, decompose=decompose, outfile='mmpbsa.in', startFrame=1, endFrame=nframe, interval=1, temperature=300, igbValue=2, name='Calculate', indi=indi, exdi=exdi)
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
        """
        The function is used to run the gmx_MMPBSA.py script
        
        Args:
          verbose: verbose level. Defaults to 0
        """
        #print("="*80)
        logging.info('Run the gmx_MMPBSA: %s'%self.mode)
        cmd = '{gmx_MMPBSA} -i {mmpbsa} -cs {complexfile} -ci {indexfile} -ct {trajectoryfile} -cp {topolfile} -cg {receptor} {ligand} -nogui >mmpbsa.log 2>&1 '.format(**self.paras)
        RC = os.system(cmd)
        if RC != 0:
            raise Exception('ERROR run: %s \nPlease ckeck the log file for details: %s'%(cmd, os.path.abspath("mmpbsa.log")))
        shutil.copy('FINAL_RESULTS_MMPBSA.dat', self.cwd)
        outfile = "FINAL_RESULTS_MMPBSA.dat"
        self.clean(verbose=verbose)
        #print('='*80)
        print('Results: %s'%outfile)

    def clean(self, verbose=0):
        """
        Clean the results
        
        Args:
          verbose: if True, print out the command line before executing it. Defaults to 0
        """
        #print('='*80)
        logging.info('Clean the results.')
        if not verbose:
            cmd = '{gmx_MMPBSA} --clean >>mmpbsa.log 2>&1 '.format(gmx_MMPBSA=gmx_MMPBSA)
            RC = os.system(cmd)
            shutil.rmtree(self.workdir)
        os.chdir(self.cwd)

    def extract_result(self, energyfile='FINAL_RESULTS_MMPBSA.dat'):
        """
        Extract the GB and PB energy from the final results file
        
        Args:
          energyfile: the name of the file that contains the energy information. Defaults to
        FINAL_RESULTS_MMPBSA.dat
        
        Returns:
          A dictionary with the GB and PB delta G values.
        """
        tagName = False#["GENERALIZED BORN", "POISSON BOLTZMANN"]
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
    parser.add_argument('-m', dest='MODE', help='Method to calculate: gb, pb, pb+gb. default:gb', default='gb', choices=['gb', 'pb', 'pb+gb', 'gb+pb'])
    parser.add_argument('-t', dest='TRAJ', help='A trajectory file contains many structure frames. File format: xtc, pdb, gro...', default=None)
    parser.add_argument('-indi', help='External dielectric constant. detault: 1.0', default=1.0, type=float)
    parser.add_argument('-dec', dest='dec', help='Decompose the energy. default:false', action='store_true', default=False)
    parser.add_argument('-D', dest='DEBUG', help='DEBUG model, keep all the files.', default=False, action='store_true')

    args = parser.parse_args()
    #exit(0)
    complexFile, topolFile, indexFile, trajFile, debug, dec, mode, indi = args.INP, args.TOP, args.ndx, args.TRAJ, args.DEBUG, args.dec, args.MODE, args.indi
    if trajFile is None:
        trajFile = complexFile
    pbsa = PBSA(mode=mode)
    pbsa.set_paras(complexfile=complexFile, trajectoryfile=trajFile, topolfile=topolFile, indexfile=indexFile, decompose=dec, indi=indi)
    pbsa.run(verbose=debug)
    detal_G = pbsa.extract_result()
    print("mode    detal_G(kcal/mole)    Std. Dev.")
    for k, v in detal_G.items():
        print('%4s    %18.4f    %9.4f'%(k, v[0], v[1]))


if __name__ == "__main__":
    main()