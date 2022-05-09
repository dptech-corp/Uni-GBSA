import os
import sys

import shutil
import parmed as pmd

from hmtpbsa.settings import GMXEXE, OMP_NUM_THREADS
from hmtpbsa.simulation.mdrun import GMXEngine
from hmtpbsa.simulation.utils import convert_to_mol2, guess_filetype, write_position_restrain

def build_lignad(ligandfile, forcefield="gaff2", charge_method="bcc", engine="acpype", clean=True, outtop=None, outcoord=None):
    """
    Build a ligand topology and coordinate file from a ligand file using acpype
    
    Args:
      ligandfile: The name of the ligand file.
      forcefield: The forcefield to use. Currently only "gaff2" is supported. Defaults to gaff2
      charge_method: The charge method to use. Valid options are "bcc" (the default) and "cm2". Defaults
    to bcc
      engine: acpype, openbabel. Defaults to acpype
      clean: If True, the temporary directory will be deleted after the calculation is finished.
    Defaults to True
    
    Returns:
      The acpype.log file is being returned.
    """
    ligandfile = os.path.abspath(ligandfile)
    ligandName = os.path.split(ligandfile)[-1][:-4]+'.TOP'
    filetype = guess_filetype(ligandfile)
    acceptFileTypes = ('pdb', 'mol2', 'mol')
    if filetype not in acceptFileTypes:
        ligandfile = convert_to_mol2(ligandfile, filetype)
    if not os.path.exists(ligandName): 
        os.mkdir(ligandName)
    cwd = os.getcwd()
    os.chdir(ligandName)
    paras = {
        'thread':OMP_NUM_THREADS,
        'ligandfile': ligandfile,
        'forcefield': forcefield,
        'method': charge_method
    }
    cmd = "export OMP_NUM_THREADS={thread};acpype -i {ligandfile} -b MOL -a {forcefield} -c {method} -f >acpype.log 2>&1 ".format(**paras)
    RC = os.system(cmd)
    if RC != 0:
        print(cmd)
        os.system('tail -n 50 acpype.log')
        raise Exception('ERROR run the acpype. see the %s for details.'%os.path.abspath("acpype.log"))
    os.chdir('MOL.acpype')
    moltop = pmd.load_file('MOL_GMX.top')
    molgro = pmd.load_file('MOL_GMX.gro', structure=True)
    if outtop:
        moltop.write(outtop)
    os.chdir(cwd)
    if clean:
        shutil.rmtree(ligandName)
    return moltop, molgro

def build_protein(pdbfile, forcefield='amber99sb-ildn', outtop=None, outcoord=None):
    """
    Build a protein topology and coordinate file from a PDB file
    
    Args:
      pdbfile: The name of the PDB file to be used as input.
      forcefield: The forcefield to use. Options are 'amber03', 'amber94', 'amber96', 'amber99',
    'amber99sb', 'amber99sb-ildn', 'amber99sb-star', 'amber14sb'. Defaults to amber99sb-ildn
    
    Returns:
      a tuple of two objects: a topology object and a Structure object.
    """
    #forcefield = {1:"amber03", 2:"amber94", 3:"amber96", 4:"amber99", 5:"amber99sb", 6:"amber99sb-ildn",
    #        7:"amber99sb-star-ildn-mut", 8:"amber14sb"}
    proteinName = os.path.split(pdbfile)[-1][:-4]+'.TOP'
    if not os.path.exists(proteinName):
        os.mkdir(proteinName)
    pdbfile = os.path.abspath(pdbfile)
    cwd = os.getcwd()
    os.chdir(proteinName)
    paras = {
        'gmx':GMXEXE,
        'pdbfile':pdbfile,
        'outfile': '1-pdb2gmx.pdb',
        'topolfile': 'topol.top',
        'forcefield': forcefield,
    }
    cmd = '{gmx} pdb2gmx -f {pdbfile} -o {outfile} -p {topolfile} -ff {forcefield} -ignh > gromacs.log 2>&1 <<EOF\n1\nEOF'.format(**paras)
    RC = os.system(cmd)
    if RC != 0:
        print(cmd)
        os.system('tail -n 50 gromacs.log')
        raise Exception('ERROR run gmx! see the log file for details %s'%os.path.abspath("gromacs.log"))
    
    engine = GMXEngine()
    boxpdb = engine.gmx_box('1-pdb2gmx.pdb', boxtype='triclinic', boxsize=0.9)
    solpdb = engine.gmx_solvate(boxpdb, 'topol.top', maxsol=5)
    ionspdb = engine.gmx_ions(solpdb, 'topol.top', conc=None, nNA=1, nCL=1, neutral=False)
    #engine.
    protgro = pmd.load_file('1-pdb2gmx.pdb', structure=True)
    #load_position_restraints('topol.top')
    prottop  = pmd.load_file('topol.top')
    if outtop:
        prottop.write(outtop)
    if outcoord:
        if outcoord[-4:] != '.pdb':
            outcoord = outcoord[:-4] + '.pdb'
        protgro.write_pdb(outcoord)

    os.chdir(cwd)
    shutil.rmtree(proteinName)
    return prottop, protgro

def build_topol(receptor, ligand, outpdb, outtop, proteinforce='amber99sb-ildn', ligandforce='gaff2'):
    """
    Build a topology file for a protein-ligand system
    
    Args:
      receptor: the name of the receptor file
      ligand: the ligand file name
      outpdb: the output pdb file name
      outtop: the output topology file name
      proteinforce: the forcefield to use for the protein (default: amber99sb-ildn). Defaults to
    amber99sb-ildn
      ligandforce: The force field to use for the ligand. Default is GAFF2. Defaults to gaff2
    """
    if isinstance(receptor, str):
        prottop, protgro = build_protein(receptor, forcefield=proteinforce)
    else:
        prottop, protgro = receptor

    if isinstance(ligand, str):
        moltop, molgro = build_lignad(ligand, forcefield=ligandforce)
    elif ligand:
        moltop, molgro = ligand

    if ligand is None:
        systop, sysgro = prottop, protgro
    else:
        systop = moltop + prottop
        sysgro = molgro + protgro
    systop.write(outtop)
    sysgro.write_pdb(outpdb)
    lines = []
    CF = 0
    records = ('NA', 'CL', 'SOL')
    with open(outtop) as fr:
        for line in fr:
            if line.startswith('[') and line.split()[1]=='molecules':
                CF = 3
                lines.append('\n; Include Position restraint file\n#ifdef POSRES\n\n#endif\n')
            if line.startswith(records) and CF:
                tmp = line.strip().split()
                molname, molnum = tmp[0], int(tmp[1])-1
                if molname == 'SOL':
                    molnum -= 2
                if molnum > 1:
                    line = '%s        %d\n'%(molname, molnum)
                else:
                    line = ''
                CF -= 1
            lines.append(line)
    with open(outtop, 'w') as fw:
        for line in lines:
            fw.write(line)
    write_position_restrain(outtop)

def main():
    pdbfile, ligandfile = sys.argv[1], sys.argv[2]

    grofile = 'sys.pdb'
    topfile = 'sys.top'

    build_topol(pdbfile, ligandfile, outpdb=grofile, outtop=topfile)

    engine = GMXEngine()
    engine.run_to_md(grofile, topfile)

if __name__  == "__main__":
    main()
