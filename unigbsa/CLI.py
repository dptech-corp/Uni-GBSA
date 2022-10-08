import os
import shutil
import sys
import argparse
from glob import glob

from .utils import generate_index_file, process_pbc

from .simulation import topology, mdrun
from .settings import logging
from .version import __version__



def PBC_main():
    parser = argparse.ArgumentParser(description='Auto process PBC for gromacs MD trajector.')
    parser.add_argument('-s', dest='tpr', help='TPR file generated from gromacs or coordinate file.', required=True)
    parser.add_argument('-f', dest='xtc', help='Trajector file to process PBC.', required=True)
    parser.add_argument('-o', dest='out', help='Result file after processed PBC.', default=None)
    parser.add_argument('-n', dest='ndx', help='Index file contains the center and output group.', default=None)
    parser.add_argument('-v', '--version', action='version', version="{prog}s ({version})".format(prog="%(prog)", version=__version__))

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
    parser.add_argument('-nt', dest='thread', help='Number of thread to run this simulation.', default=4)
    parser.add_argument('-verbose', help='Keep the directory or not.', default=False, action='store_true')
    parser.add_argument('-v', '--version', action='version', version="{prog}s ({version})".format(prog="%(prog)", version=__version__))

    args = parser.parse_args()
    protein, ligand, outdir, cF = args.protein, args.ligand, args.outdir, args.c
    proteinForcefield, ligandForcefield = args.protforce, args.ligforce
    verbose = args.verbose
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
        outtop, outcoord, outitp = os.path.join(outdir, ligandName+'.top'), os.path.join(outdir, ligandName+'.pdb'), os.path.join(outdir, ligandName+'.itp')
        # outcoord parameter is useless
        ligtop, liggro = topology.build_lignad(ligandfile, forcefield=ligandForcefield, charge_method='bcc', outtop=outtop, outcoord=outcoord, itpfile=outitp, nt=args.thread)
        if cF:
            comxtop, comxcoord = os.path.join(outdir, "%s_%s.top"%(proteinName, ligandName)), os.path.join(outdir, "%s_%s.pdb"%(proteinName, ligandName))
            topology.build_topol((prottop, protgro), (ligtop, liggro), outtop=comxtop, outpdb=comxcoord, verbose=verbose, nt=args.thread)

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
    parser.add_argument('-nt', dest='thread', help='Number of thread to run this simulation.', default=4)
    parser.add_argument('-v', '--version', action='version', version="{prog}s ({version})".format(prog="%(prog)", version=__version__))

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
        topology.build_topol(receptor, None, outpdb=grofile, outtop=topfile, nt=args.thread)

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
            topology.build_topol(receptor, ligandfile, outpdb=grofile, outtop=topfile, ligandforce=ligandForcefield, nt=args.thread)

            logging.info('Building simulation for: %s'%ligandName)
            engine = mdrun.GMXEngine()
    
            ionspdb, topfile = engine.run_to_ions(grofile, topfile, rundir=None, boxtype=boxtype, boxsize=boxsize, conc=conc)

            shutil.copy(ionspdb, "complex.pdb")
            shutil.copy(topfile, "complex.top")

            engine.clean(pdbfile=grofile)
    os.chdir(cwd)

def simulation_run():
    parser = argparse.ArgumentParser(description='Run MD simulation for input file.')
    parser.add_argument('-p', dest='protein', help='Protein file for the simulation.', required=True)
    parser.add_argument('-l', dest='ligand', help='Ligand file or directory for the simulation.', default="")
    parser.add_argument('-pf', dest='protforce', help='Protein forcefield.', default='amber03')
    parser.add_argument('-lf', dest='ligforce', help='Ligand forcefiled: gaff or gaff2.', default='gaff', choices=['gaff','gaff2'])
    parser.add_argument('-bt', dest='boxtype', help='Simulation box type, default: triclinic', default='triclinic')
    parser.add_argument('-box', help='Simulation box size.', nargs=3, type=float, default=None)
    parser.add_argument('-d', help='Distance between the solute and the box.', default=0.9, type=float)
    parser.add_argument('-conc', help='Specify salt concentration (mol/liter). default=0.15', default=0.15, type=float)
    parser.add_argument('-o', dest='outdir', help='A output directory.', default=None)
    parser.add_argument('-nsteps', dest='nstep', help='Simulation steps. default:2500', default=2500, type=int)
    parser.add_argument('-nframe', dest='nframe', help='Number of frame to save for the xtc file. default:100', default=100, type=int)
    parser.add_argument('-nt', dest='thread', help='Number of thread to run this simulation.', default=4)
    parser.add_argument('-verbose', help='Keep all the files in the simulation.', action='store_true', default=False)
    parser.add_argument('-v', '--version', action='version', version="{prog}s ({version})".format(prog="%(prog)", version=__version__))    

    args = parser.parse_args()
    proteinfile, ligand, outdir = args.protein, args.ligand, args.outdir
    proteinForcefield, ligandForcefield = args.protforce, args.ligforce
    boxtype, box, conc, boxsize, nsteps, nframe, nt = args.boxtype, args.box, args.conc, args.d, args.nstep, args.nframe, args.thread
    verbose = args.verbose
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
        topology.build_topol(receptor, None, outpdb=grofile, outtop=topfile, nt=nt)

        logging.info('Build simulation for %s'%proteinName)
        engine = mdrun.GMXEngine()
        mdgro, mdxtc, topfile = engine.run_to_md(grofile, topfile, rundir=None, boxtype=boxtype, boxsize=boxsize, conc=conc, nsteps=nsteps, nframe=nframe, nt=nt)

        shutil.copy(mdgro, os.path.join(outdir, '%s_system.gro'%proteinName))
        shutil.copy(topfile, os.path.join(outdir, '%s_system.top'%proteinName))
        shutil.copy(mdxtc, os.path.join(outdir, '%s_traj.xtc'%proteinName))
        if not verbose:
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
            topology.build_topol(receptor, ligandfile, outpdb=grofile, outtop=topfile, ligandforce=ligandForcefield, nt=nt)

            logging.info('Building simulation for: %s'%ligandName)
            engine = mdrun.GMXEngine()
    
            mdgro, mdxtc, topfile = engine.run_to_md(grofile, topfile, rundir=None, boxtype=boxtype, boxsize=boxsize, conc=conc, nsteps=nsteps, nframe=nframe, nt=nt)

            shutil.copy(mdgro, "complex.gro")
            shutil.copy(topfile, "complex.top")
            shutil.copy(mdxtc, 'traj_comx.xtc')
            if not verbose:
                engine.clean(pdbfile=grofile)
    os.chdir(cwd)

def traj_pipeline(args=None):
    from unigbsa.gbsa.pbsarun import PBSA
    parser = argparse.ArgumentParser(description='Free energy calcaulated by PBSA method.')
    parser.add_argument('-i', dest='INP', help='A pdb file or a tpr file for the trajectory.', required=True)
    parser.add_argument('-p', dest='TOP', help='Gromacs topol file for the system.', required=True)
    parser.add_argument('-ndx', dest='ndx', help='Gromacs index file, must contain recepror and ligand group.', required=True)
    parser.add_argument('-m', dest='mode', help='MMPBSA mode', nargs='+', default=['GB'])
    parser.add_argument('-f', dest='mmpbsafile', help='Input mmpbsa file', default=None)
    parser.add_argument('-t', dest='TRAJ', help='A trajectory file contains many structure frames. File format: xtc, pdb, gro...', default=None)
    parser.add_argument('-nt', dest='thread', help='Set number of thread to run this program.', type=int, default=1)
    parser.add_argument('-D', dest='DEBUG', help='DEBUG model, keep all the files.', default=False, action='store_true')
    parser.add_argument('-v', '--version', action='version', version="{prog}s ({version})".format(prog="%(prog)", version=__version__))
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)
    #exit(0)
    complexFile, topolFile, indexFile, trajFile, debug, mmpbsafile, nt = args.INP, args.TOP, args.ndx, args.TRAJ, args.DEBUG, args.mmpbsafile, args.thread
    if trajFile is None:
        trajFile = complexFile
    if mmpbsafile:
        mmpbsafile = os.path.abspath(mmpbsafile)
        pbsaParas = None
    else:
        pbsaParas = { "modes":','.join(args.mode)}

    pbsa = PBSA()
    pbsa.set_paras(complexfile=complexFile, trajectoryfile=trajFile, topolfile=topolFile, indexfile=indexFile, mmpbsafile=mmpbsafile, pbsaParas=pbsaParas, nt=nt)
    pbsa.run(verbose=debug)
    detal_G = pbsa.extract_result()
    print("mode    detal_G(kcal/mole)    Std. Dev.")
    for k, v in detal_G.items():
        print('%4s    %18.4f    %9.4f'%(k, v[0], v[1]))

def mmpbsa_plot():
    from unigbsa.gbsa import plots
    parser = argparse.ArgumentParser(description='Analysis and plot results for MMPBSA.')
    parser.add_argument('-i', help='MMPBSA result directory. Which contains FINAL_RESULTS_MMPBSA.dat, FINAL_DECOMP_MMPBSA.dat, EO.csv or DEO.csv file.', required=True)
    parser.add_argument('-o', help='Figure output directory. default: analysis', default='analysis')
    parser.add_argument('-v', '--version', action='version', version="{prog}s ({version})".format(prog="%(prog)", version=__version__))    


    args = parser.parse_args()
    inp, oup = args.i, args.o
    filenames = ["FINAL_RESULTS_MMPBSA.dat", "FINAL_DECOMP_MMPBSA.dat", "EO.csv", "DEO.csv"]
    funcdic = {
        "FINAL_RESULTS_MMPBSA.dat": plots.analysis_FINAL,
        "FINAL_DECOMP_MMPBSA.dat": plots.analysis_DECOMP,
        "EO.csv": plots.analysis_traj_EO,
        "DEO.csv": plots.analysis_traj_DEO
    }
    for fname in filenames:
        datfile = os.path.join(inp, fname)
        if os.path.exists(datfile):
            func = funcdic[fname]
            logging.info('Analysis file: %s'%fname)
            func(datfile, outdir=oup)