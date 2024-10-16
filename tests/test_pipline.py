import os
import shutil
import unittest
from unigbsa.pipeline import md_pipeline, minim_pipeline, traj_pipeline, load_configue_file

import warnings


TEST_EM_CONFIG = os.path.dirname(os.path.abspath(__file__)) + '/config/em.ini'
TEST_MD_CONFIG = os.path.dirname(os.path.abspath(__file__)) + '/config/md.ini'


class TestPipline(unittest.TestCase):
    def setUpClass():
        warnings.simplefilter('ignore', ResourceWarning)

    def base(self, pdbfile, ligandfile):
        pdbfile = os.path.abspath(pdbfile)
        ligandfile = os.path.abspath(ligandfile)
        pdbname = "/tmp/" + os.path.split(pdbfile)[-1]
        if os.path.exists(pdbname):
            shutil.rmtree(pdbname)
        os.mkdir(pdbname)
        return pdbfile, ligandfile, pdbname

#################### version 0.0.2 ########################

    def pipeline_simulation(self, pdbfile, ligandfiles, configfile=None, nt=1):
        pdbfile, ligandfile, workdir = self.base(pdbfile, ligandfiles[0])
        cwd = os.getcwd()
        os.chdir(workdir)
        paras = load_configue_file(configfile)
        if paras['simulation']['mode']=='em':
            minim_pipeline(receptorfile=pdbfile, ligandfiles=[ligandfile], paras=paras, outfile='BindingEnergy.csv', nt=nt, verbose=True)
        elif paras['simulation']['mode']=='md':
            md_pipeline(receptorfile=pdbfile, ligandfiles=[ligandfile], paras=paras, outfile='BindingEnergy.csv', nt=nt)
        EF = os.path.exists('BindingEnergy.csv')
        self.assertTrue(EF)
        os.chdir(cwd)
        shutil.rmtree(workdir)

    def pipeline_minima(self, pdbfile, ligandfiles):
        self.pipeline_simulation(pdbfile, ligandfiles, TEST_EM_CONFIG)
    
    def pipeline_md(self, pdbfile, ligandfiles):
        self.pipeline_simulation(pdbfile, ligandfiles, TEST_MD_CONFIG)

    def test_traj(self, nt=1):
        pdbfile, topfile, indexfile = '../example/3f/complex.pdb', os.path.abspath('../example/3f/complex.top'), os.path.abspath('../example/3f/index.ndx')
        pdbfile, ligandfile, workdir = self.base(pdbfile, pdbfile)
        cwd = os.getcwd()
        os.chdir(workdir)
        paras = load_configue_file()
        traj_pipeline(pdbfile, pdbfile, topfile, indexfile, pbsaParas=paras['GBSA'])
        EF = os.path.exists('Energy.csv')
        self.assertTrue(EF)
        os.chdir(cwd)
        shutil.rmtree(workdir)

    def test_minima(self):
        pdbfile = '../example/1ceb/1ceb_protein.pdb'
        ligandfiles = ['../example/1ceb/1ceb_ligand.sdf']
        self.pipeline_minima(pdbfile, ligandfiles)

    def test_md(self):
        pdbfile = '../example/1ceb/1ceb_protein.pdb'
        ligandfiles = ['../example/1ceb/1ceb_ligand.sdf']
        self.pipeline_md(pdbfile, ligandfiles)

    def test_protein_protein(self):
        pdbfile = '../example/protein-protein/1vf6_REC.pdb'
        ligandfiles = ['../example/protein-protein/1vf6_LIG.pdb']
        self.pipeline_minima(pdbfile, ligandfiles)

if __name__ == "__main__":
    unittest.main()
