import os
import sys
import shutil
import logging
import argparse

import pandas as pd


from pbsa import PBSA
from utils import generate_index_file
from mdrun.topbuild import build_topol, GMXEngine, build_protein


LOG_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
DATE_FORMAT = "%m/%d/%Y %H:%M:%S %p"

logging.basicConfig(level=logging.INFO, format=LOG_FORMAT, datefmt=DATE_FORMAT)

def traj_pipeline(complexfile, trajfile, topolfile, indexfile, mode='gb', dec=False, debug=False):
    pbsa = PBSA(mode=mode)
    pbsa.set_paras(complexfile=complexfile, trajectoryfile=trajfile, topolfile=topolfile, indexfile=indexfile, decompose=dec)
    pbsa.run(verbose=debug)
    detal_G = pbsa.extract_result()
    print("mode    detal_G(kcal/mole)    Std. Dev.")
    for k, v in detal_G.items():
        print('%4s    %18.4f    %9.4f'%(k, v[0], v[1]))
    return detal_G

def base_pipeline(receptorfile, ligandfile, paras=None):
    grofile = 'complex.pdb'
    topfile = 'complex.top'
    build_topol(receptorfile, ligandfile, outpdb=grofile, outtop=topfile)
    pass

def minim_peipline(receptorfile, ligandfiles, paras, outfile='BindingEnergy.csv'):
    simParas = paras['simulation']
    gbsaParas = paras['GBSA']

    receptorfile = os.path.abspath(receptorfile)
    logging.info('Build protein topology.')
    receptor = build_protein(receptorfile, forcefield=simParas['proteinforcefield'])

    detalGdict = {'name':[]}
    for k in gbsaParas['mode'].split('+'):
        detalGdict[k.upper()] = []
        detalGdict['%s_err'%k.upper()] = []
    
    cwd = os.getcwd()
    for ligandfile in ligandfiles:
        ligandfile = os.path.abspath(ligandfile)
        ligandName = os.path.split(ligandfile)[-1][:-4]
        if not os.path.exists(ligandName):
            os.mkdir(ligandName)
        os.chdir(ligandName)

        grofile = 'complex.pdb'
        topfile = 'complex.top'
        logging.info('Build ligand topology: %s'%ligandName)
        build_topol(receptor, ligandfile, outpdb=grofile, outtop=topfile, ligandforce=simParas['ligandforcefield'])

        logging.info('Running energy minimization: %s'%ligandName)
        engine = GMXEngine()

        minimgro, outtop = engine.run_to_minim(grofile, topfile, boxtype=simParas['boxtype'], boxsize=simParas['boxsize'], conc=simParas['conc'])
    
        cmd = 'gmx editconf -f %s -o %s >/dev/null 2>&1'%(minimgro, grofile)
        RC = os.system(cmd)
        if RC!=0:
            raise Exception('Error conver %s to %s'%(minimgro, grofile))
        shutil.copy(topfile, outtop)

        indexfile = generate_index_file(grofile)
        gbsaParas = paras['GBSA']
        detalG = traj_pipeline(grofile, trajfile=grofile, topolfile=topfile, indexfile=indexfile, mode=gbsaParas['mode'], dec=gbsaParas['decomposition'])
        detalGdict['name'].append(ligandName)
        for k,v in detalG.items():
            detalGdict[k.upper()].append(v[0])
            detalGdict["%s_err"%k.upper()].append(v[1])
        engine.clean(pdbfile=grofile)
        os.chdir(cwd)
    df = pd.DataFrame(detalGdict)
    df.to_csv(outfile, index=False)

def md_pipeline(receptorfile, ligandfiles, paras, outfile='BindingEnergy.csv'):
    simParas = paras['simulation']
    gbsaParas = paras['GBSA']

    receptorfile = os.path.abspath(receptorfile)
    logging.info('Build protein topology.')
    receptor = build_protein(receptorfile, forcefield=simParas['proteinforcefield'])
    detalGdict = {'name':[]}
    for k in gbsaParas['mode'].split('+'):
        detalGdict[k.upper()] = []
        detalGdict['%s_err'%k.upper()] = []

    cwd = os.getcwd()
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
        build_topol(receptor, ligandfile, outpdb=grofile, outtop=topfile, ligandforce=simParas['ligandforcefield'])

        logging.info('Running simulation: %s'%ligandName)
        engine = GMXEngine()
    
        mdgro, mdxtc, outtop = engine.run_to_md(grofile, topfile, boxtype=simParas['boxtype'], boxsize=simParas['boxsize'], conc=simParas['conc'], nstep=simParas['nstep'])

        cmd = 'gmx editconf -f %s -o %s >/dev/null 2>&1'%(mdgro, grofile)
        RC = os.system(cmd)
        if RC!=0:
            raise Exception('Error conver %s to %s'%(mdgro, grofile))

        shutil.copy(topfile, outtop)
        shutil.copy(mdxtc, xtcfile)

        logging.info('Running GBSA: %s'%ligandName)
        indexfile = generate_index_file(grofile)
    
        detalG = traj_pipeline(grofile, trajfile=xtcfile, topolfile=topfile, indexfile=indexfile, mode=gbsaParas['mode'], dec=gbsaParas['decomposition'])
        detalGdict['name'].append(ligandName)
        for k,v in detalG.items():
            detalGdict[k.upper()].append(v[0])
            detalGdict["%s_err"%k.upper()].append(v[1])
        engine.clean(pdbfile=grofile)
        os.chdir(cwd)
    df = pd.DataFrame(detalGdict)
    df.to_csv(outfile, index=False)


def main():
    parser = argparse.ArgumentParser(description='GBSA Calculation.')
    parser.add_argument('-i', dest='receptor', help='Input protein file with pdb format.', required=True)
    parser.add_argument('-l', dest='ligand', help='Ligand files to calculate binding energy.', nargs='+', required=True)
    parser.add_argument('-c', dest='config', help='Configue file', default='')

    args = parser.parse_args()
    receptor, ligands, conf = args.receptor, args.ligand, args.config
    #base_pipeline(receptorfile=receptor, ligandfile=ligands[0], paras=conf)
    paras = {
        'simulation':{
            'boxtype' : 'triclinic',
            'boxsize' :0.9,
            'conc': 0.15,
            'nstep':500000,
            'proteinforcefield': 'amber99sb-ildn',
            'ligandforcefield': 'gaff2'
        },
        'GBSA':{
            'mode': 'gb',
            'decomposition': False
        }
    }
    #minim_peipline(receptorfile=receptor, ligandfiles=ligands, paras=paras)
    md_pipeline(receptorfile=receptor, ligandfiles=ligands, paras=paras)

if __name__ == "__main__":
    main()