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

#################### version 0.0.1 ########################

    def pipline_single(self, pdbfile, rID, lID):
        pdbname = os.path.split(pdbfile)[-1]
        #cmd = 'hmtpbsa -i %s -o /tmp/%s > hmtpbsa.log 2>&1 << EOF \n %d,%d \n EOF'%(pdbfile, pdbname, rID, lID)
        cmd = 'hmtpbsa -i %s -o /tmp/%s  << EOF \n %d,%d \n EOF'%(pdbfile, pdbname, rID, lID)
        print(cmd)
        os.system(cmd)
        EF = os.path.exists('/tmp/%s/FINAL_RESULTS_MMPBSA.dat'%pdbname)
        self.assertTrue(EF)
        os.system('rm -rf /tmp/%s/'%pdbname)
    
    def pipline_traj(self, pdbfile, trajfile, rID, lID):
        pdbname = os.path.split(pdbfile)[-1]
        #cmd = 'hmtpbsa -i %s -o /tmp/%s > hmtpbsa.log 2>&1 << EOF \n %d,%d \n EOF'%(pdbfile, pdbname, rID, lID)
        cmd = 'hmtpbsa -i %s -t %s -o /tmp/%s  << EOF \n %d,%d \n EOF'%(pdbfile, trajfile, pdbname, rID, lID)
        print(cmd)
        os.system(cmd)
        EF = os.path.exists('/tmp/%s/FINAL_RESULTS_MMPBSA.dat'%pdbname)
        self.assertTrue(EF)
        os.system('rm -rf /tmp/%s/'%pdbname)

    #def test_2fvy(self):
    #    self.pipline_single('../example/2fvy.pdb', 1, 2)
    
    #def test_traj_ST(self):
    #    self.pipline_traj('../example/ST/com.tpr', '../example/ST/com_traj.xtc', 1, 2)

#################### version 0.0.2 ########################

    def pipeline_simulation(self, pdbfile, ligandfiles, configfile):
        pdbfile, ligandfile, workdir = self.base(pdbfile, ligandfiles[0])
        ligand = ' '.join(ligandfiles)
        cwd = os.getcwd()
        os.chdir(workdir)
        cmd = 'export NUM_OF_THREAD=1;hmtpbsa-pipeline -i %s -l %s -c %s'%(pdbfile, ligandfile, configfile)
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
        pdbfile = '../example/md/protein.pdb'
        ligandfiles = ['../example/md/3f.mol']
        self.pipeline_minima(pdbfile, ligandfiles)

    def test_md(self):
        pdbfile = '../example/md/protein.pdb'
        ligandfiles = ['../example/md/3f.mol']
        self.pipeline_md(pdbfile, ligandfiles)
if __name__ == "__main__":
    unittest.main()
