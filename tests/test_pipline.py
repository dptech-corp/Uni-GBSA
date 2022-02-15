import os
import unittest

class TestPipline(unittest.TestCase):
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

    def test_2fvy(self):
        self.pipline_single('../example/2fvy.pdb', 1, 2)
    
    def test_traj_ST(self):
        self.pipline_traj('../example/ST/com.tpr', '../example/ST/com_traj.xtc', 1, 2)

if __name__ == "__main__":
    unittest.main()
