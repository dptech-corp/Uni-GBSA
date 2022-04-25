import os
import shutil
import sys
import argparse
from glob import glob

from .utils import generate_index_file, process_pbc

from .simulation import topology, mdrun
from .settings import logging

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

def topol_builder():
    parser = argparse.ArgumentParser(description='Build topology file for input file.')
    parser.add_argument('-p', dest='protein', help='Protein file or directory to build topology.', default="")
    parser.add_argument('-l', dest='ligand', help='Ligand file or directory to build topology.', default="")
    parser.add_argument('-pf', dest='protforce', help='Protein forcefield.', default='amber03')
    parser.add_argument('-lf', dest='ligforce', help='Ligand forcefiled: gaff or gaff2.', default='gaff', choices=['gaff','gaff2'])
    parser.add_argument('-o', dest='outdir', help='A output directory.', default='GMXtop')
    parser.add_argument('-c', help='Combine the protein and ligand topology. Suppport for one protein and more ligands. default:True', action='store_true', default=True)
    
    args = parser.parse_args()
    protein, ligand, outdir, cF = args.protein, args.ligand, args.outdir, args.c
    proteinForcefield, ligandForcefield = args.protforce, args.ligforce
    if not protein and not ligand:
        print('Not found input file!')

    proteinfiles, ligandfiles = [], []
    if os.path.isdir(protein):
        proteinfiles = glob(os.path.join(protein, '*.pdb'))
    elif os.path.isfile(protein):
        proteinfiles = [protein]
    
    if os.path.isdir(ligand):
        ligandfiles = glob(os.path.join(ligand, "*"))
    elif os.path.isfile(ligand):
        ligandfiles = [ligand]
    
    if len(proteinfiles) >=2 or len(ligandfiles)==0 or len(proteinfiles)==0:
        cF = False
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outdir = os.path.abspath(outdir)

    for proteinfile in proteinfiles:
        proteinName = os.path.split(proteinfile)[-1][:-4]
        outtop, outcoord = os.path.join(outdir, proteinName+'.top'), os.path.join(outdir, proteinName+'.pdb')
        prottop, protgro = topology.build_protein(proteinfile, forcefield=proteinForcefield, outtop=outtop, outcoord=outcoord)
    for ligandfile in ligandfiles:
        ligandName = os.path.split(ligandfile)[-1][:-4]
        outtop, outcoord = os.path.join(outdir, ligandName+'.top'), os.path.join(outdir, ligandName+'.pdb')
        # outcoord parameter is useless
        ligtop, liggro = topology.build_lignad(ligandfile, forcefield=ligandForcefield, charge_method='bcc', outtop=outtop, outcoord=outcoord)
        if cF:
            comxtop, comxcoord = os.path.join(outdir, "%s_%s.top"%(proteinName, ligandName)), os.path.join(outdir, "%s_%s.pdb"%(proteinName, ligandName))
            topology.build_topol((prottop, protgro), (ligtop, liggro), outtop=comxtop, outpdb=comxcoord)

def simulation_builder():
    parser = argparse.ArgumentParser(description='Build MD simulation for input file.')
    parser.add_argument('-p', dest='protein', help='Protein file for the simulation.', required=True)
    parser.add_argument('-l', dest='ligand', help='Ligand file or directory for the simulation.', default="")
    parser.add_argument('-pf', dest='protforce', help='Protein forcefield.', default='amber03')
    parser.add_argument('-lf', dest='ligforce', help='Ligand forcefiled: gaff or gaff2.', default='gaff', choices=['gaff','gaff2'])
    parser.add_argument('-bt', dest='boxtype', help='Simulation box type, default: triclinic', default='triclinic')
    parser.add_argument('-box', help='Simulation box size.', nargs=3, type=float, default=None)
    parser.add_argument('-d', help='Distance between the solute and the box.', default=0.9, type=float)
    parser.add_argument('-conc', help='Specify salt concentration (mol/liter). default=0.15', default=0.15, type=float)
    parser.add_argument('-o', dest='outdir', help='A output directory.', default=None)

    args = parser.parse_args()
    proteinfile, ligand, outdir = args.protein, args.ligand, args.outdir
    proteinForcefield, ligandForcefield = args.protforce, args.ligforce
    boxtype, box, conc, boxsize = args.boxtype, args.box, args.conc, args.d
    if box:
        boxsize = box

    ligandfiles = []
    if os.path.isdir(ligand):
        ligandfiles = glob(os.path.join(os.path.abspath(ligand), "*"))
    elif os.path.isfile(ligand):
        ligandfiles = [ os.path.abspath(ligand) ]

    proteinName = os.path.split(proteinfile)[-1][:-4]
    receptor = topology.build_protein(proteinfile, forcefield=proteinForcefield)

    if not outdir:
        outdir = 'GMXtop'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outdir = os.path.abspath(outdir)

    cwd = os.getcwd()
    os.chdir(outdir)
    grofile = 'complex.pdb'
    topfile = 'complex.top'

    if len(ligandfiles) == 0:
        logging.info('No ligand found, build protein only.')
        topology.build_topol(receptor, None, outpdb=grofile, outtop=topfile)

        logging.info('Build simulation for %s'%proteinName)
        engine = mdrun.GMXEngine()
        ionspdb, topfile = engine.run_to_ions(grofile, topfile, rundir=None, boxtype=boxtype, boxsize=boxsize, conc=conc)

        shutil.copy(ionspdb, os.path.join(outdir, '%s_system.pdb'%proteinName))
        shutil.copy(topfile, os.path.join(outdir, '%s_system.top'%proteinName))
        engine.clean(pdbfile=grofile)
        os.system('rm complex.pdb complex.top >/dev/null 2>&1')
    else:
        for ligandfile in ligandfiles:
            ligandfile = os.path.abspath(ligandfile)
            ligandName = os.path.split(ligandfile)[-1][:-4]
            if not os.path.exists(ligandName):
                os.mkdir(ligandName)
            os.chdir(ligandName)

            logging.info('Build ligand topology: %s'%ligandName)
            topology.build_topol(receptor, ligandfile, outpdb=grofile, outtop=topfile, ligandforce=ligandForcefield)

            logging.info('Building simulation for: %s'%ligandName)
            engine = mdrun.GMXEngine()
    
            ionspdb, topfile = engine.run_to_ions(grofile, topfile, rundir=None, boxtype=boxtype, boxsize=boxsize, conc=conc)

            shutil.copy(ionspdb, "complex.pdb")
            shutil.copy(topfile, "complex.top")

            engine.clean(pdbfile=grofile)
    os.chdir(cwd)