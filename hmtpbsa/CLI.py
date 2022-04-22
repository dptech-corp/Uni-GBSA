import os
import sys
import argparse
from .utils import generate_index_file, process_pbc

def PBC_main():
    parser = argparse.ArgumentParser(description='Auto process PBC for gromacs MD trajector.')
    parser.add_argument('-s', dest='tpr', help='TPR file generated from gromacs or coordinate file.', required=True)
    parser.add_argument('-f', dest='xtc', help='Trajector file to process PBC.', required=True)
    parser.add_argument('-o', dest='out', help='Result file after processed PBC.', default=None)
    parser.add_argument('-n', dest='ndx', help='Index file contains the center and output group.', default=None)

    args = parser.parse_args()
    tprfile, trajfile, outfile, indexfile = args.tpr, args.xtc, args.out, args.ndx
    if outfile is None:
        suffix = trajfile[-4:]
        outfile = os.path.spplit(trajfile)[-1][:-4] + '-pbc' + suffix
    if indexfile is None:
        indexfile = generate_index_file(tprfile, pbc=True)
    process_pbc(trajfile, tprfile, indexfile=indexfile, outfile=outfile)
