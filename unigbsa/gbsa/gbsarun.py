import os
import shutil
import pandas as pd

from unigbsa.settings import gmx_MMPBSA, MPI, logging
from unigbsa.utils import obtain_id_from_index
from unigbsa.gbsa.utils import obtain_num_of_frame
from unigbsa.gbsa.parameters import generate_input_file
from unigbsa.gbsa.io import parse_GMXMMPBSA_RESULTS


'''
1. one protein file and many other ligands file
2. one complex file and a trajectory file.
'''

class GBSA(object):
    def __init__(self, workdir='UNIGBSA') -> None:
        self.workdir = os.path.abspath(workdir)
        self.cwd = os.getcwd()
        self.verbose = 0
        self.deltaG = None
        self.resG = None

    def set_paras(self, complexfile, trajectoryfile, topolfile, indexfile, pbsaParas=None, mmpbsafile=None, nt=1):
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
        self.nframe = obtain_num_of_frame(trajectoryfile)
        if trajectoryfile == complexfile:
            ext = trajectoryfile[-3:]
            os.system('cp %s com_traj.%s'%(trajectoryfile, ext))
            trajectoryfile = os.path.abspath('com_traj.%s'%ext)
        if mmpbsafile is None:
            mmpbsafile = generate_input_file(pbsaParas, outfile='mmpbsa.in')
            #shutil.copy(mmpbsafile, '../mmpbsa.in')
        mmpbsafile = os.path.abspath(mmpbsafile)
            
        # mode='gb', outfile='mmpbsa.in', startFrame=1, endFrame=1, interval=1, temperature=300, igbValue=2, name='Calculate', decompose=False
        receptor, ligand = obtain_id_from_index(indexfile)
        self.paras = {
            'mmpbsa': mmpbsafile,
            'gmx_MMPBSA':gmx_MMPBSA,
            'complexfile': complexfile,
            'indexfile': indexfile,
            'trajectoryfile': trajectoryfile,
            'topolfile': topolfile,
            'receptor': receptor,
            'ligand': ligand,
            'eofile':'FINAL_EO_MMPBSA.csv',
            'deofile':'FINAL_DEO_MMPBSA.csv',
            'numthread':nt
        }
        return mmpbsafile

    def run(self, verbose=0):
        """
        The function is used to run the gmx_MMPBSA.py script
        
        Args:
          verbose: verbose level. Defaults to 0
        """
        #print("="*80)
        logging.info('Run the MMPB(GB)SA.')
        if self.nframe > 2 and MPI:
            cmd = 'mpirun --use-hwthread-cpus --allow-run-as-root -np {numthread} \
                {gmx_MMPBSA} MPI -i {mmpbsa} -cs {complexfile} -ci {indexfile} -ct {trajectoryfile} -cp {topolfile} -cg {receptor} {ligand} -nogui >mmpbsa.log 2>&1 '.format(**self.paras)
        else:
            cmd = '{gmx_MMPBSA} MPI -i {mmpbsa} -cs {complexfile} -ci {indexfile} -ct {trajectoryfile} -cp {topolfile} -cg {receptor} {ligand} -nogui >mmpbsa.log 2>&1 '.format(**self.paras)
        RC = os.system(cmd)
        if RC != 0:
            raise Exception('ERROR run: %s \nPlease ckeck the log file for details: %s'%(cmd, os.path.abspath("mmpbsa.log")))
        self.save_results()
        self.verbose = verbose
        self.clean(verbose=verbose)
        print('='*80)
        print('Results: Energy.csv Dec.csv')

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

    def extract_result(self, energyfile='Energy.csv'):
        if self.deltaG is None:
             self.save_results()
        return self.deltaG

    def extract_result_v14(self, energyfile='FINAL_RESULTS_MMPBSA.dat'):
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

    def extract_result_v15(self, energyfile='FINAL_RESULTS_MMPBSA.dat'):
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
                elif line.startswith('ΔTOTAL'):
                    if tagName:
                        lineTemp = line.split()
                        detal_G[tagName] = (float(lineTemp[1]), float(lineTemp[2]))
                    else:
                        logging.warning("Found a DELTA G without name!")
        return detal_G
        
    def save_results(self, mmxsafile=None):
        if mmxsafile is None:
            mmxsafile = os.path.join(self.workdir, 'COMPACT_MMXSA_RESULTS.mmxsa')
        if not os.path.exists(mmxsafile):
            logging.warning('Not found mmxsa file!')
            return
        self.deltaG, self.resG = parse_GMXMMPBSA_RESULTS(mmxsafile=mmxsafile)
