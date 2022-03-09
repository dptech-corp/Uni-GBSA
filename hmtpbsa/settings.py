import os
import logging

LOG_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
DATE_FORMAT = "%m/%d/%Y %H:%M:%S %p"
DEFAULT_CONFIGURE_FILE = os.path.dirname(os.path.abspath(__file__))+ '/data/detault.ini'

# base configure
GMXEXE='gmx'
MDPFILESDIR = os.path.dirname(os.path.abspath(__file__)) + '/simulation/mdp'


logging.basicConfig(level=logging.INFO, format=LOG_FORMAT, datefmt=DATE_FORMAT)


gmx_MMPBSA='gmx_MMPBSA'
if 'AMBERHOME' not in os.environ:
    #os.environ['AMBERHOME'] = set_amber_home(gmx_MMPBSA)
    raise Exception("Not found variable AMBERHOME")

if 'NUM_OF_THREAD' in os.environ:
    NUM_OF_THREAD = os.environ['NUM_OF_THREAD']
else:
    NUM_OF_THREAD = 4