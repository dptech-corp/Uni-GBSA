import os

import json
import shutil
import argparse
import itertools
import pandas as pd
import multiprocessing

from copy import deepcopy as copy
from concurrent.futures import ProcessPoolExecutor

from unigbsa.simulation.mdrun import GMXEngine
from unigbsa.settings import PathManager, GMXEXE
from unigbsa.simulation.topology import build_topol, build_protein
from unigbsa.pipeline import traj_pipeline
from unigbsa.utils import load_configue_file, generate_index_file
from unigbsa.utils import logging


KEY = ['ligandName', 'Frames', 'mode', 'complex','receptor','ligand','Internal','Van der Waals','Electrostatic','Polar Solvation','Non-Polar Solvation','Gas','Solvation','TOTAL']


def reres_gro(infile, outfile):
    cmd = '%s editconf -f %s -o %s -resnr 1 >/dev/null 2>&1'%(GMXEXE, infile, outfile)
    RC = os.system(cmd)
    if RC!=0:
        raise Exception('Error convert %s to %s'%(infile, outfile))
    return outfile

def threads_split(njob, nt):
    threads = max([round(nt/njob), 1])
    nt = njob if njob < nt else nt
    return threads, nt

def load_scan_paras(jsonfile: str, scantype='fixed') -> dict:
    parasdict = {}
    with open(jsonfile)  as fr:
        data = json.load(fr)
    defaultparas = load_configue_file()
    scantypes = 'fixed all'
    simulationparas = {
        'forcefield':{
            "proteinforcefield": "amber03",
            "ligandforcefield": "gaff",
            "ligandCharge": "bcc"
        }, 
        'simulation':{
            "mode": "em",
            "boxtype": "triclinic",
            "boxsize": 0.9,
            "conc": 0.15,
            "nsteps": 500000,
            "nframe": 100,
            "eqsteps": 50000,
            "maxsol": 0,
        }
    }

    modedict = {
        'md': 'nsteps',
        'gb': 'igb',
        'pb': 'ipb'
    }
    varparas = {}
    for ki, vi in data.items():
        for kj, vi in vi.items():
            if not isinstance(vi, list):
                defaultparas[ki][kj] = vi
            elif isinstance(vi, list) and len(vi) == 1:
                defaultparas[ki][kj] = vi[0]
            elif (ki=="simulation" and kj in ['mode', 'ligandCharge', 'ligandforcefield', 'proteinforcefield']) or \
                 (ki == "GBSA" and kj in ['modes', 'indi', 'exdi']):
                key = f'{ki}_{kj}'
                varparas[key] = vi
    if scantype == 'fixed':
        for k, v in varparas.items():
            ki, kj = k.split('_')
            for vi in v:
                dic = copy(defaultparas)
                if 'mode' in kj and '-' in vi:
                    mode, modevalue = vi.split('-')
                    dic[ki][kj] = mode
                    if mode == 'md':
                        modevalue = int(modevalue)
                    dic[ki][modedict[mode]] = modevalue
                else:
                    dic[ki][kj] = vi
                name = '%s_%s'%(str(k), str(vi))
                parasdict[name] = dic
    elif scantype == 'all':
        keys = varparas.keys()
        values = [varparas[k] for k in keys]
        groups = itertools.product(*values)
        for group in groups:
            dic = copy(defaultparas)
            name = []
            for k,v in zip(keys, group):
                ki, kj = k.split('_')
                if 'mode' in kj and '-' in v:
                    mode, modevalue = v.split('-')
                    dic[ki][kj] = mode
                    if mode == 'md':
                        modevalue = int(modevalue)
                    dic[ki][modedict[mode]] = modevalue
                else:
                    dic[ki][kj] = v
                name.append('%s_%s'%(str(k), str(v)))
            parasdict['__'.join(name)] = dic
    else:
        raise Exception(f'scantype {scantype} not one of the: {scantypes}')

    parasdict_unique = {}
    for k, v in parasdict.items():
        if v not in parasdict_unique.values():
            parasdict_unique[k] = v

    simulationparas = {}
    for name, v in parasdict_unique.items():
        simparas = v['simulation']
        forcefieldkey = f"{simparas['proteinforcefield']}_{simparas['ligandforcefield']}_{simparas['ligandCharge']}"
        if simparas['mode'] == 'md':
            simulationkey = f"{simparas['mode']}_{simparas['nsteps']}"
        else:
            simulationkey = f"{simparas['mode']}"
        if forcefieldkey not in simulationparas:
            simulationparas[forcefieldkey] = {simulationkey:v}
        elif simulationkey not in simulationparas[forcefieldkey]:
            simulationparas[forcefieldkey][simulationkey] = v

    return parasdict_unique, simulationparas

def calc_R2(expfile, gbsa):
    exp = pd.read_csv(expfile)
    if isinstance(gbsa, str):
        GBSA = pd.read_csv(gbsa)
    else:
        GBSA = gbsa
    df = pd.merge(exp, GBSA, on='ligandName')
    R  = df[['TOTAL', 'dG_exp']].corr().loc['TOTAL', 'dG_exp']
    return R, R**2

def iter_paras(args) -> None:
    '''
    '''
    from unigbsa.pipeline import minim_pipeline, md_pipeline, base_pipeline
    receptor, ligands, name, paras, outfile, expdatfile, verbose, nt = args
    
    cwd = os.getcwd()
    if not os.path.exists(name):
        os.mkdir(name)
    os.chdir(name)
    if paras['simulation']['mode'] == 'em':
        minim_pipeline(receptorfile=receptor, ligandfiles=ligands, paras=paras, outfile=outfile, verbose=verbose, nt=nt)
    elif paras['simulation']['mode'] == 'md':
        md_pipeline(receptorfile=receptor, ligandfiles=ligands, paras=paras, outfile=outfile, verbose=verbose, nt=nt)
    elif paras['simulation']['mode'] == 'input':
        base_pipeline(receptorfile=receptor, ligandfiles=ligands, paras=paras, outfile=outfile, verbose=verbose, nt=nt)
    GBSA = pd.read_csv(outfile)
    with open(name+'.json', 'w') as fw:
        json.dump(paras, fw, indent=4)
    os.chdir(cwd)
    outfile = os.path.join(name, outfile)
    R, R2 = calc_R2(expdatfile, outfile)
    
    return (name, R, R2)

def build_topology_walker(arg):
    receptor, ligandfile, paras, nt = arg
    paras = copy(paras)
    simParas = paras['simulation']
    ligandName = os.path.split(ligandfile)[-1][:-4]
    with PathManager(ligandName) as pm:
        grofile = 'complex.pdb'
        topfile = 'complex.top'
        logging.info('Build ligand topology: %s'%ligandName)
        try:
            build_topol(receptor, ligandfile, outpdb=grofile, outtop=topfile, proteinforce=simParas['proteinforcefield'], ligandforce=simParas['ligandforcefield'], charge_method=simParas['ligandCharge'], nt=nt)
            grofile = pm.abspath(grofile)
            topfile = pm.abspath(topfile)
        except Exception as e:
            print(e)
            logging.warning('Failed to generate forcefield for ligand: %s'%ligandName)
            return None
    outfiles = {  
            'complexfile': grofile,
            'topolfile0': topfile,
            'GBSAinput': grofile
        }
    return ligandName, outfiles

def build_topology_MPI(receptorfiles, ligandfiles, paras, outdir, nt=4):
    threads, nworker = threads_split(len(ligandfiles), nt)
    with PathManager(outdir) as pm:
        if isinstance(receptorfiles, str):
            receptor = build_protein(pm.abspath(receptorfiles, parent=True), forcefield=paras['simulation']['proteinforcefield'])
            receptors = [receptor] * len(ligandfiles)
        else:
            receptors = pm.abspath(receptorfiles, parent=True)
        ligandfiles = pm.abspath(ligandfiles, parent=True)
        args = [ (receptor, ligandfile, paras, threads) for receptor, ligandfile in zip(receptors, ligandfiles) ]
        with ProcessPoolExecutor(max_workers=nworker) as pool:
            outdict = { out[0]:out[1] for out in list(pool.map(build_topology_walker, args)) if out is not None }
    outparas = copy(paras)
    outparas['files'] = outdict
    return outparas

def structural_optimization_walker(arg):
    paras, ligandName, outdir, threads = arg
    paras = copy(paras)
    files = paras['files'][ligandName]
    simParas = paras['simulation']
    complexfile = files['complexfile']
    topolfile = files['topolfile0']
    if outdir is None:
        ligandir = os.path.split(complexfile)[0]
    else:
        ligandir = os.path.join(outdir, ligandName)
    with PathManager(ligandir) as pm:
        engine = GMXEngine()
        if paras['simulation']['mode'] == 'em':
            GBSAInputfile = 'complex_minim.pdb'
            minimgro, outtop = engine.run_to_minim(complexfile, topolfile, boxtype=simParas['boxtype'], boxsize=simParas['boxsize'], conc=simParas['conc'], maxsol=simParas['maxsol'], nt=threads)
            GBSAInputfile = reres_gro(minimgro, 'complex_minim.pdb')
            files['GBSAinput'] = pm.abspath(GBSAInputfile)
            files['GBSAtraj'] = pm.abspath(GBSAInputfile)
            shutil.copy(outtop, 'topol.top')
            files['topolfile'] = pm.abspath('topol.top')
            engine.clean(pdbfile=complexfile)
        elif paras['simulation']['mode'] == 'md':
            mdgro, mdxtc, outtop = engine.run_to_md(complexfile, topolfile, boxtype=simParas['boxtype'], boxsize=simParas['boxsize'], conc=simParas['conc'], nsteps=simParas['nsteps'], nframe=simParas['nframe'], eqsteps=simParas['eqsteps'], nt=threads)
            GBSAInputfile = reres_gro(mdgro, 'complex_md.pdb')
            files['GBSAinput'] = pm.abspath(GBSAInputfile)
            files['GBSAtraj'] = pm.abspath('traj_comx.xtc')
            shutil.copy(outtop, 'topol.top')
            files['topolfile'] = pm.abspath('topol.top')
            shutil.copy(mdxtc, 'traj_comx.xtc')
            engine.clean(pdbfile=complexfile)
        elif paras['simulation']['mode'] == 'input':
            files['GBSAinput'] = complexfile
            files['GBSAtraj'] = complexfile
            files['topolfile'] = topolfile
        else:
            return None
        indexfile = generate_index_file(files['GBSAinput'])
        files['indexfile'] = indexfile
    return ligandName, files

def structural_optimization_MPI(paras, outdir=None, nt=4):
    ligandNames = paras['files'].keys()
    threads, nworker = threads_split(len(ligandNames), nt)
    args = [ (paras, ligandName, outdir, threads) for ligandName in ligandNames ]
    with ProcessPoolExecutor(max_workers=nworker) as pool:
        outfiles = { out[0]:out[1] for out in list(pool.map(structural_optimization_walker, args)) if out is not None }
    outparas = copy(paras)
    outparas['files'] = outfiles
    return outparas

def gbsa_calculation_walker(arg):
    paras, ligandName, outdir, threads = arg
    files = paras['files'][ligandName]
    complexfile = files['complexfile']
    if outdir is None:
        ligandir = os.path.split(complexfile)[0]
    else:
        ligandir = os.path.join(outdir, ligandName)
    with PathManager(ligandir) as pm:
        deltaG = traj_pipeline(files['GBSAinput'], trajfile=files['GBSAtraj'], topolfile=files['topolfile'], indexfile=files['indexfile'], pbsaParas=paras['GBSA'], mmpbsafile=None, nt=threads, verbose=False)
        deltaG['ligandName'] = ligandName
    return ligandName, deltaG

def gbsa_calculation_MPI(paras, outdir, nt=4):
    ligandNames = paras['files'].keys()
    threads, nworker = threads_split(len(ligandNames), nt)
    args = [ (paras, ligandName, outdir, threads) for ligandName in ligandNames ]
    with ProcessPoolExecutor(max_workers=nworker) as pool:
        results = [ out for out in list(pool.map(gbsa_calculation_walker, args)) if out is not None ]
    df = None
    for result in results:
        if df is None:
            df = result[1]
        else:
            df = pd.concat([df, result[1]])
    outparas = copy(paras)
    if df is None:
        logging.warn('Faild run gbsa for paras: %s'%str(paras))
        return None
    df = df.reset_index(drop=True)
    outparas['results'] = df.to_dict()
    parafile = os.path.join(outdir, 'paras'+'.json')
    with open(parafile, 'w') as fw:
        json.dump(outparas, fw, indent=4)
    
    outparas['results'] = df
    outcsv = os.path.join(outdir, 'deltaG.csv')
    df.to_csv(outcsv, index=False)
    return outparas

def scan_parameters_v2(receptors, protdir, ligands, ligdir, expdatfile, parasfile, outdir, nt=4) -> None:
    '''
    '''
    if ligands is None:
        ligands = []
    if receptors is None:
        receptors = []
    if not isinstance(receptors, list):
        receptors = [receptors]
    if protdir:
        for fileName in os.listdir(protdir):
            if fileName.endswith('.pdb'):
                receptors.append(os.path.join(protdir, fileName))
    if ligdir:
        for fileName in os.listdir(ligdir):
            if fileName.endswith(('mol', 'sdf')):
                ligands.append(os.path.join(ligdir, fileName))
    if len(ligands) == 0 or len(receptors) == 0:
        raise Exception('No ligands or receptors file found.')
    if len(receptors) == 1 and len(ligands) != 1:
        receptors = receptors * len(ligands)
    with PathManager(outdir) as pm:
        expdatfile = pm.abspath(expdatfile, parent=True)
        parasfile = pm.abspath(parasfile, parent=True)
        ligands = pm.abspath(sorted([lig for lig in ligands]), parent=True)
        receptors = pm.abspath(sorted([prot for prot in receptors]), parent=True)
        logging.info('load scan paras.')
        parasdicts, simulationparas = load_scan_paras(parasfile)
        for name, parasdic in simulationparas.items():
            logging.info('Building protein and ligand topology.')
            topfileparas = build_topology_MPI(receptors, ligands, parasdic[list(parasdic.keys())[0]], name, nt=nt)
            for k, v in parasdic.items():
                outdir = os.path.join(name, k)
                topfileparas['simulation'] = v['simulation']
                outparas = structural_optimization_MPI(topfileparas, outdir=outdir, nt=nt)
                simulationparas[name][k] = outparas
        outset = []
        for k, v in parasdicts.items():
            simparas = v['simulation']
            forcefieldkey = f"{simparas['proteinforcefield']}_{simparas['ligandforcefield']}_{simparas['ligandCharge']}"
            if simparas['mode'] == 'md':
                simulationkey = f"{simparas['mode']}_{simparas['nsteps']}"
            else:
                simulationkey = f"{simparas['mode']}"
            paras = simulationparas[forcefieldkey][simulationkey]
            paras['GBSA'] = v['GBSA']
            modes = paras['GBSA']['modes']
            if '-' in modes:
                gmode, gtype = modes.split('-')
                paras['GBSA']['modes'] = gmode
                if gmode.upper() == 'GB':
                    paras['GBSA']['igb'] = gtype
                elif gmode.upper() == 'PB':
                    paras['GBSA']['ipb'] = gtype
            outdir = os.path.join(forcefieldkey, f"{simulationkey}/{k}")
            logging.info(f'GBSA calculation: {forcefieldkey} {simulationkey} {k}')
            outparas = gbsa_calculation_MPI(paras, outdir=outdir, nt=nt)
            if outparas is None:
                R, R2 = 0, -10
            else:
                R, R2 = calc_R2(expdatfile, outparas['results'])
            logging.info(f'GBSA results: {R} {R2}')
            outparas['results']['R'] = R
            outparas['results']['R2'] = R2
            outset.append((k, R, R2, pm.abspath(outdir+'/paras.json')))
        logging.info('Write output: paras_performance.csv')
        df = pd.DataFrame(outset, columns=['name', 'R', 'R2', 'parasjson'])
        df = df.sort_values(by='R2', ascending=False).reset_index(drop=True)
        R2max = list(df.loc[0, :])
        df.to_csv('paras_performance.csv', index=False)
    print('='*80)
    print('The best para name is: %s'%R2max[0])
    print('The best para R2 is: %.4f'%R2max[2])
    print('The best para file is: %s'%R2max[3])
    print('='*80)

def scan_parameters(receptor, ligands, ligdir, expdatfile, parasfile, verbose, outdir, nt=4) -> None:
    '''
    '''
    outfile = 'deltaG.csv'
    if ligands is None:
        ligands = []
    if ligdir:
        for fileName in os.listdir(ligdir):
            if fileName.endswith(('mol', 'sdf')):
                ligands.append(os.path.join(ligdir, fileName))
    if len(ligands)==0:
        raise Exception('No ligands file found.')
    ligands = [ os.path.abspath(l) for l in ligands ]
    receptor = os.path.abspath(receptor)

    if not os.path.exists(outdir):
        os.mkdir(outdir)
    cwd = os.getcwd()
    os.chdir(outdir)

    R2max = ('', 0)
    parasdict = load_scan_paras(parasfile)
    args = [ (receptor, ligands, k, v, outfile, expdatfile, verbose, 1) for k,v in parasdict.items() ]
    with ProcessPoolExecutor(max_workers=nt) as pool:
        outset = list(pool.map(iter_paras, args))
    with open('paras_performance.csv', 'w') as fw:
        fw.write('name,R,R2\n')
        for outline in outset:
            fw.write('%s,%.6f,%.6f\n'%(outline[0], outline[1], outline[2]))
            if outline[2]>R2max[1]:
                R2max = (outline[0], outline[2])
    os.chdir(cwd)
    print('The best para name is: %s'%R2max[0])
    print('The best para R2 is: %.4f'%R2max[1])

class ParameterScan(object):
    threads = 1
    paras = None
    def __init__(self) -> None:
        pass

def main():
    parser = argparse.ArgumentParser(description='Perform an automatic parameter optimization prior to production MM/GB(PB)SA calculations.')
    parser.add_argument('-i', dest='receptor', help='Input protein file in pdb format.', default=None)
    parser.add_argument('-pd', dest='protdir', help='Directory containing many protein files. file format: .pdb', default=None)
    parser.add_argument('-l', dest='ligand', help='Ligand files to calculate binding energy.', nargs='+', default=None)
    parser.add_argument('-ld', dest='ligdir', help='Directory containing many ligand files. file format: .mol or .sdf', default=None)
    parser.add_argument('-e', help='Experiment data file.', required=True)
    parser.add_argument('-c', dest='parasfile', help='Parameters to scan', required=True)
    parser.add_argument('-o', dest='outdir', help='Output directory.', default='pbsa.scan')
    parser.add_argument('-nt', dest='threads', help='Set number of threads to run this program.', type=int, default=multiprocessing.cpu_count())
    parser.add_argument('--verbose', help='Keep all the files.', action='store_true', default=False)

    args = parser.parse_args()
    receptor, protdir, ligands, ligdir, expdatfile, parasfile, verbose, outdir, nt = args.receptor, args.protdir, args.ligand, args.ligdir, os.path.abspath(args.e), os.path.abspath(args.parasfile), args.verbose, args.outdir, args.threads
    # scan_parameters(receptor, ligands, ligdir, expdatfile, parasfile, verbose, outdir, nt)
    #if isinstance(receptor, str):
    #    receptor = [receptor] * len(ligands)
    scan_parameters_v2(receptor, protdir, ligands, ligdir, expdatfile, parasfile, outdir, nt=nt)
    
