import os
import re
import sys
import uuid
import shutil
import logging

LOG_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
DATE_FORMAT = "%m/%d/%Y %H:%M:%S %p"
DEFAULT_CONFIGURE_FILE = os.path.dirname(os.path.abspath(__file__))+ '/data/default.ini'
TEMPLATE_TOP = os.path.dirname(os.path.abspath(__file__))+ '/data/template.json'
PBSA_PARAMETER_FILE = os.path.dirname(os.path.abspath(__file__))+ '/data/mmpbsa.in'

# base configure
def find_gmx():
    RC = os.system('gmx -h >/dev/null 2>&1')
    if RC == 0:
        return 'gmx'
    RC = os.system('gmx_mpi -h >/dev/null 2>&1')
    if RC == 0:
        return 'gmx_mpi'
    logging.error('Not found gmx or gmx_mpi.')
    exit(1)


    
GMXEXE = find_gmx()
MDPFILESDIR = os.path.dirname(os.path.abspath(__file__)) + '/simulation/mdp'

def has_mpirun():
    RC = os.system('which mpirun >/dev/null 2>&1')
    if RC == 0:
        return True
    else:
        return False

MPI = has_mpirun()

logging.basicConfig(level=logging.INFO, format=LOG_FORMAT, datefmt=DATE_FORMAT)


gmx_MMPBSA='gmx_MMPBSA'
if 'AMBERHOME' not in os.environ:
    #os.environ['AMBERHOME'] = set_amber_home(gmx_MMPBSA)
    print("Not found variable AMBERHOME")

if 'OMP_NUM_THREADS' in os.environ:
    OMP_NUM_THREADS = os.environ['OMP_NUM_THREADS']
else:
    OMP_NUM_THREADS = 1

def set_OMP_NUM_THREADS(nt):
    os.environ['OMP_NUM_THREADS'] = str(nt)

def obtain_MMPBSA_version():
    versionFile = '/tmp/' + uuid.uuid1().hex
    os.system('gmx_MMPBSA -v >%s 2>&1 '%versionFile)
    with open(versionFile) as fr:
        text = fr.read()
        p = re.compile('v(\d\.\d)\.\d', re.S)
        m = re.search(p, text)
    os.system('rm %s'%versionFile)
    if m:
        version = float(m.groups()[0])
    else:
        logging.error('Can not found gmx_MMPBSA version.')
        print(text)
        sys.exit(1)
    return version

PBSA_VERSION = 1.5  #obtain_MMPBSA_version()

class PathManager(object):
    def __init__(self, path) -> None:
        self.path = path
        self.cwd = os.getcwd()
    
    def __enter__(self) -> str:
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        os.chdir(self.path)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        os.chdir(self.cwd)

    def abspath(self, files, parent=False):
        '''
        Obtain absolute file path by parent directory
        '''
        if parent:
            cwd = self.cwd
        else:
            cwd = '.'
        if isinstance(files, str):
            return os.path.abspath(os.path.join(cwd, files))
        else:
            return [os.path.abspath(os.path.join(cwd, f)) for f in files]
