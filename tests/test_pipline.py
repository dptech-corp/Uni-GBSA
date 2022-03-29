import os
import shutil
import unittest

TEST_EM_CONFIG = os.path.dirname(os.path.abspath(__file__)) + '/config/em.ini'
TEST_MD_CONFIG = os.path.dirname(os.path.abspath(__file__)) + '/config/md.ini'
class TestPipline(unittest.TestCase):
    def base(self, pdbfile, ligandfile):
        pdbfile = os.path.abspath(pdbfile)
        ligandfile = os.path.abspath(ligandfile)
        pdbname = "/tmp/" + os.path.split(pdbfile)[-1]
        if os.path.exists(pdbname):
            shutil.rmtree(pdbname)
        os.mkdir(pdbname)
        return pdbfile, ligandfile, pdbname

#################### version 0.0.2 ########################

    def pipeline_simulation(self, pdbfile, ligandfiles, configfile):
        pdbfile, ligandfile, workdir = self.base(pdbfile, ligandfiles[0])
        ligand = ' '.join(ligandfiles)
        cwd = os.getcwd()
        os.chdir(workdir)
        cmd = 'export OMP_NUM_THREADS=1;hmtpbsa-pipeline -i %s -l %s -c %s'%(pdbfile, ligandfile, configfile)
        print(cmd)
        os.system(cmd)
        EF = os.path.exists('BindingEnergy.csv')
        self.assertTrue(EF)
        os.chdir(cwd)
        shutil.rmtree(workdir)

    def pipeline_minima(self, pdbfile, ligandfiles):
        self.pipeline_simulation(pdbfile, ligandfiles, TEST_EM_CONFIG)
    
    def pipeline_md(self, pdbfile, ligandfiles):
        self.pipeline_simulation(pdbfile, ligandfiles, TEST_MD_CONFIG)

    def test_traj(self):
        pdbfile, topfile, indexfile = '../example/3f/complex.pdb', os.path.abspath('../example/3f/complex.top'), os.path.abspath('../example/3f/index.ndx')
        pdbfile, ligandfile, workdir = self.base(pdbfile, pdbfile)
        cwd = os.getcwd()
        os.chdir(workdir)
        cmd = 'hmtpbsa-traj -i %s -p %s -ndx %s -m pb+gb -t %s -indi 1'%(pdbfile, topfile, indexfile, pdbfile)
        RC = os.system(cmd)
        if RC != 0:
            EF = False
        else:
            EF = True
        self.assertTrue(EF)
        os.chdir(cwd)
        shutil.rmtree(workdir)

    def test_minima(self):
        pdbfile = '../example/2fvy/protein.pdb'
        ligandfiles = ['../example/2fvy/BGC.mol2']
        self.pipeline_minima(pdbfile, ligandfiles)

    def test_md(self):
        pdbfile = '../example/2fvy/protein.pdb'
        ligandfiles = ['../example/2fvy/BGC.mol2']
        self.pipeline_md(pdbfile, ligandfiles)

if __name__ == "__main__":
    unittest.main()
