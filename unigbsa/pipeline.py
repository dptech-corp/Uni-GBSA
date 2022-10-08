import os
import shutil
import argparse
import traceback
import pandas as pd

from unigbsa.version import __version__
from unigbsa.gbsa.gbsarun import GBSA
from unigbsa.utils import generate_index_file, load_configue_file
from unigbsa.simulation.mdrun import GMXEngine
from unigbsa.simulation.topology import build_topol, build_protein
from unigbsa.settings import logging, DEFAULT_CONFIGURE_FILE, GMXEXE, set_OMP_NUM_THREADS

KEY = ['ligandName', 'Frames', 'mode', 'complex','receptor','ligand','Internal','Van der Waals','Electrostatic','Polar Solvation','Non-Polar Solvation','Gas','Solvation','TOTAL', 'status']
def traj_pipeline(complexfile, trajfile, topolfile, indexfile, pbsaParas=None, mmpbsafile=None, nt=1, verbose=False):
    """
    A pipeline for calculate GBSA/PBSA for trajectory
    
    Args:
      complexfile: The name of the PDB file containing the protein-ligand complex.
      trajfile: the trajectory file, in pdb format
      topolfile: The topology file of the complex.
      indexfile: the index file for the complex
      mode: gb or tm. Defaults to gb
      dec: whether to decompose the free energy into its components. Defaults to False
      debug: if True, will print out all the debug messages. Defaults to False
    
    Returns:
      detal_G is a dictionary, the key is the mode, the value is a list, the first element is the
    average value, the second element is the standard deviation.
    """

    reresfile = complexfile[:-4]+'_reres.pdb'
    cmd = '%s editconf -f %s -o %s -resnr 1 >/dev/null 2>&1'%(GMXEXE, complexfile, reresfile)
    RC = os.system(cmd)
    if RC!=0:
       raise Exception('Error convert %s to %s'%(complexfile, reresfile))
    pbsa = GBSA()
    mmpbsafile = pbsa.set_paras(complexfile=reresfile, trajectoryfile=trajfile, topolfile=topolfile, indexfile=indexfile, pbsaParas=pbsaParas, mmpbsafile=mmpbsafile, nt=nt)
    pbsa.run(verbose=verbose)
    detal_G = pbsa.extract_result()
    print("Frames    mode    detal_G(kcal/mole)")
    for i, irow in detal_G.iterrows():
        print('%6d    %4s    %18.4f  '%(irow['Frames'], irow['mode'], irow['TOTAL']))
    return detal_G

def base_pipeline(receptorfile, ligandfiles, paras, nt=1, mmpbsafile=None, outfile='BindingEnergy.csv', verbose=False):
    """
    This function takes a receptorfile and ligandfile, and build a complex.pdb and complex.top file
    
    Args:
      receptorfile: the file name of the receptor pdb file
      ligandfile: the name of the ligand file
      paras: a dictionary of parameters for the pipeline.
    """
    simParas = paras['simulation']
    pbsaParas = paras['GBSA']

    receptorfile = os.path.abspath(receptorfile)
    logging.info('Build protein topology.')
    receptor = build_protein(receptorfile, forcefield=simParas['proteinforcefield'])
    
    cwd = os.getcwd()
    df = None
    ligandnames = []
    status = []
    pbsaParas = paras['GBSA']
    d = pd.DataFrame({'Frames': 1, 'mode':pbsaParas['modes'], 'complex':0.0,'receptor':0.0,'ligand':0.0,'Internal':0.0,'Van der Waals':0.0,'Electrostatic':0,'Polar Solvation':0.0,'Non-Polar Solvation':0.0,'Gas':0.0,'Solvation':0.0,'TOTAL':0.0}, index=[1])
    for ligandfile in ligandfiles:
        statu = 'S'
        ligandfile = os.path.abspath(ligandfile)
        ligandName = os.path.split(ligandfile)[-1][:-4]
        if not os.path.exists(ligandName):
            os.mkdir(ligandName)
        os.chdir(ligandName)

        grofile = 'complex.pdb'
        topfile = 'complex.top'
        logging.info('Build ligand topology: %s'%ligandName)
        try:
            build_topol(receptor, ligandfile, outpdb=grofile, outtop=topfile, ligandforce=simParas['ligandforcefield'], charge_method=simParas['ligandCharge'], nt=nt)
        except Exception as e:
            if len(ligandfiles)==1:
                traceback.print_exc()
            statu = 'F_top'
            dl = d
            logging.warning('Failed to generate forcefield for ligand: %s'%ligandName)

        indexfile = generate_index_file(grofile)
        
        if statu == 'S':
            try:
                dl = traj_pipeline(grofile, trajfile=grofile, topolfile=topfile, indexfile=indexfile, pbsaParas=pbsaParas, mmpbsafile=mmpbsafile, verbose=verbose, nt=nt)
            except:
                if len(ligandfiles)==1:
                    traceback.print_exc()
                statu = 'F_GBSA'
                dl = d
                logging.warning('Failed to run GBSA for ligand: %s'%ligandName)
        ligandnames.extend([ligandName]*len(dl))
        status.extend([statu]*len(dl))
        if df is None:
            df = dl
        else:
            df = pd.concat([df, dl])
        os.chdir(cwd)
    df['ligandName'] = ligandnames
    df['status'] = status
    df[KEY].to_csv(outfile, index=False)



def minim_pipeline(receptorfile, ligandfiles, paras, mmpbsafile=None, nt=1, outfile='BindingEnergy.csv', verbose=False):
    """
    It runs the simulation pipeline for each ligand.
    
    Args:
      receptorfile: The name of the receptor file.
      ligandfiles: a list of ligand files
      paras: a dictionary of parameters
      outfile: the output file name. Defaults to BindingEnergy.csv
    """
    simParas = paras['simulation']
    pbsaParas = paras['GBSA']

    receptorfile = os.path.abspath(receptorfile)
    logging.info('Build protein topology.')
    receptor = build_protein(receptorfile, forcefield=simParas['proteinforcefield'])
    
    cwd = os.getcwd()
    df = None
    ligandnames = []
    status = []
    d = pd.DataFrame({'Frames': 1, 'mode':pbsaParas['modes'], 'complex':0.0,'receptor':0.0,'ligand':0.0,'Internal':0.0,'Van der Waals':0.0,'Electrostatic':0,'Polar Solvation':0.0,'Non-Polar Solvation':0.0,'Gas':0.0,'Solvation':0.0,'TOTAL':0.0}, index=[1])
    for ligandfile in ligandfiles:
        statu = 'S'
        ligandfile = os.path.abspath(ligandfile)
        ligandName = os.path.split(ligandfile)[-1][:-4]
        if not os.path.exists(ligandName):
            os.mkdir(ligandName)
        os.chdir(ligandName)

        grofile = 'complex.pdb'
        topfile = 'complex.top'
        logging.info('Build ligand topology: %s'%ligandName)
        try:
            build_topol(receptor, ligandfile, outpdb=grofile, outtop=topfile, ligandforce=simParas['ligandforcefield'], charge_method=simParas['ligandCharge'], nt=nt)
        except Exception as e:
            if len(ligandfiles)==1:
                traceback.print_exc()
                exit(256)
            statu = 'F_top'
            dl = d
            logging.warning('Failed to generate forcefield for ligand: %s'%ligandName)   

        logging.info('Running energy minimization: %s'%ligandName)
        engine = GMXEngine()
        try:
            minimgro, outtop = engine.run_to_minim(grofile, topfile, boxtype=simParas['boxtype'], boxsize=simParas['boxsize'], conc=simParas['conc'], maxsol=simParas['maxsol'], nt=nt)
            cmd = '%s editconf -f %s -o %s -resnr 1 >/dev/null 2>&1'%(GMXEXE, minimgro, grofile)
            RC = os.system(cmd)
            if RC!=0:
                raise Exception('Error convert %s to %s'%(minimgro, grofile))
            shutil.copy(topfile, outtop)
        except Exception as e:
            if len(ligandfiles)==1:
                traceback.print_exc()
                exit(256)
            statu = 'F_md'
            dl = d
            logging.warning('Failed to run simulation for ligand: %s'%ligandName)  


        indexfile = generate_index_file(grofile)
        if statu == 'S':
            try:
                dl = traj_pipeline(grofile, trajfile=grofile, topolfile=topfile, indexfile=indexfile, pbsaParas=pbsaParas, mmpbsafile=mmpbsafile, verbose=verbose, nt=nt)
            except:
                if len(ligandfiles)==1:
                    traceback.print_exc()
                statu = 'F_GBSA'
                dl = d
                logging.warning('Failed to run GBSA for ligand: %s'%ligandName)
        ligandnames.extend([ligandName]*len(dl))
        status.extend([statu]*len(dl))
        if df is None:
            df = dl
        else:
            df = pd.concat([df, dl])
        if not verbose:
            engine.clean(pdbfile=grofile)
        os.chdir(cwd)

        os.chdir(cwd)
    df['ligandName'] = ligandnames
    df['status'] = status
    df[KEY].to_csv(outfile, index=False)

def md_pipeline(receptorfile, ligandfiles, paras, mmpbsafile=None, nt=1, outfile='BindingEnergy.csv', verbose=False):
    """
    The main function of this script
    
    Args:
      receptorfile: the protein file
      ligandfiles: a list of ligand files
      paras: a dictionary of parameters
      outfile: the output file name. Defaults to BindingEnergy.csv
    """
    simParas = paras['simulation']
    pbsaParas = paras['GBSA']

    receptorfile = os.path.abspath(receptorfile)
    logging.info('Build protein topology.')
    receptor = build_protein(receptorfile, forcefield=simParas['proteinforcefield'])

    cwd = os.getcwd()
    df = None
    ligandnames = []
    status = []
    for ligandfile in ligandfiles:
        print('='*80)
        ligandfile = os.path.abspath(ligandfile)
        ligandName = os.path.split(ligandfile)[-1][:-4]
        if not os.path.exists(ligandName):
            os.mkdir(ligandName)
        os.chdir(ligandName)

        grofile = 'complex.pdb'
        topfile = 'complex.top'
        xtcfile = 'traj_com.xtc'
        logging.info('Build ligand topology: %s'%ligandName)
        build_topol(receptor, ligandfile, outpdb=grofile, outtop=topfile, ligandforce=simParas['ligandforcefield'], charge_method=simParas['ligandCharge'], nt=nt)

        logging.info('Running simulation: %s'%ligandName)
        engine = GMXEngine()
    
        mdgro, mdxtc, outtop = engine.run_to_md(grofile, topfile, boxtype=simParas['boxtype'], boxsize=simParas['boxsize'], conc=simParas['conc'], nsteps=simParas['nsteps'], nframe=simParas['nframe'], eqsteps=simParas['eqsteps'], nt=nt)

        cmd = '%s editconf -f %s -o %s -resnr 1 >/dev/null 2>&1'%(GMXEXE, mdgro, grofile)
        RC = os.system(cmd)
        if RC!=0:
            raise Exception('Error convert %s to %s'%(mdgro, grofile))

        shutil.copy(topfile, outtop)
        shutil.copy(mdxtc, xtcfile)

        #logging.info('Running GBSA: %s'%ligandName)
        indexfile = generate_index_file(grofile)
        if 'startframe' not in pbsaParas:
            pbsaParas["startframe"] = 2
        deltaG = traj_pipeline(grofile, trajfile=xtcfile, topolfile=topfile, indexfile=indexfile, pbsaParas=pbsaParas, mmpbsafile=mmpbsafile, nt=nt, verbose=verbose)
        ligandnames.extend([ligandName]*simParas['nframe'])
        status.extend(['S']*simParas['nframe'])
        if df is None:
            df = deltaG
        else:
            df = pd.concat([df, deltaG])
        if not verbose:
            engine.clean(pdbfile=grofile)
        os.chdir(cwd)
    df['ligandName'] = ligandnames
    df['status'] = status
    df[KEY].to_csv(outfile, index=False)


def main(args=None):
    parser = argparse.ArgumentParser(description='GBSA Calculation.  Version: %s'%__version__)
    parser.add_argument('-i', dest='receptor', help='Input protein file with pdb format.', required=True)
    parser.add_argument('-l', dest='ligand', help='Ligand files to calculate binding energy.', nargs='+', default=None)
    parser.add_argument('-c', dest='config', help='Configue file, default: %s'%DEFAULT_CONFIGURE_FILE, default=DEFAULT_CONFIGURE_FILE)
    parser.add_argument('-d', dest='ligdir', help='Floder contains many ligand files. file format: .mol or .sdf', default=None)
    parser.add_argument('-f', dest='pbsafile', help='gmx_MMPBSA input file. default=None', default=None)
    parser.add_argument('-o', dest='outfile', help='Output file.', default='BindingEnergy.csv')
    parser.add_argument('-nt', dest='thread', help='Set number of thread to run this program.', type=int, default=1)
    parser.add_argument('--decomp', help='Decompose the free energy. default:False', action='store_true', default=False)
    parser.add_argument('--verbose', help='Keep all the files.', action='store_true', default=False)
    parser.add_argument('-v', '--version', action='version', version="{prog}s ({version})".format(prog="%(prog)", version=__version__))

    args = parser.parse_args(args)
    receptor, ligands, conf, ligdir, outfile, decomposition, nt, verbose = args.receptor, args.ligand, args.config, args.ligdir, args.outfile, args.decomp, args.thread, args.verbose
    
    if ligands is None:
        ligands = []
    if ligdir:
        for fileName in os.listdir(ligdir):
            if fileName.endswith(('mol','sdf')):
                ligands.append(os.path.join(ligdir, fileName))
    if len(ligands)==0:
        raise Exception('No ligands file found.')

    if not os.path.exists(conf):
        raise Exception("Not found the configure file! %s"%conf)

    mmpbsafile = os.path.abspath(args.pbsafile) if args.pbsafile else args.pbsafile
    set_OMP_NUM_THREADS(nt)
    paras = load_configue_file(conf)
    if decomposition:
        paras['GBSA']['modes'] += ',decomposition'

    if paras['simulation']['mode'] == 'em':
        minim_pipeline(receptorfile=receptor, ligandfiles=ligands, paras=paras, outfile=outfile, mmpbsafile=mmpbsafile, verbose=verbose, nt=nt)
    elif paras['simulation']['mode'] == 'md':
        md_pipeline(receptorfile=receptor, ligandfiles=ligands, paras=paras, outfile=outfile, mmpbsafile=mmpbsafile, verbose=verbose, nt=nt)
    elif paras['simulation']['mode'] == 'input':
        base_pipeline(receptorfile=receptor, ligandfiles=ligands, paras=paras, outfile=outfile, mmpbsafile=mmpbsafile, verbose=verbose, nt=nt)

if __name__ == "__main__":
    main()
