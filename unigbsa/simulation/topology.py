import os
import sys

import uuid
import json
import shutil
import parmed as pmd

from unigbsa.settings import GMXEXE, TEMPLATE_TOP
from unigbsa.simulation.mdrun import GMXEngine
from unigbsa.simulation.utils import convert_format, guess_filetype, fix_insertions
from unigbsa.simulation.utils import assign_partial_charge, write_position_restrain
from unigbsa.simulation.utils import obtain_net_charge, gen_index_for_gbsa

def build_lignad(ligandfile, forcefield="gaff2", charge_method="bcc", engine="acpype", verbose=False, outtop=None, outcoord=None, molname='MOL', itpfile=None, sqm_opt=True, nt=1):
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
    acceptFileTypes = ('mol2', 'mol')
    
    acpype_charge_methods = ['gas', 'bcc']
    if charge_method not in acpype_charge_methods:
        ligandfile = assign_partial_charge(ligandfile, filetype, charge_method=charge_method)
        charge_method = 'user'
    elif filetype != 'mol':
        ligandfile = convert_format(ligandfile, filetype)
    if not os.path.exists(ligandName):
        os.mkdir(ligandName)
    charge = obtain_net_charge(ligandfile)
    cwd = os.getcwd()
    os.chdir(ligandName)
    if sqm_opt:
        sqm_key = "grms_tol=0.005,qm_theory='AM1',scfconv=1.d-8,ndiis_attempts=700,maxcyc=0"
    else:
        sqm_key = "grms_tol=0.0005,qm_theory='AM1',scfconv=1.d-10,ndiis_attempts=700,"
    paras = {
        'thread': nt,
        'ligandfile': ligandfile,
        'forcefield': forcefield,
        'method': charge_method,
        'molname': molname,
        'net_charge': charge,
        'sqm_key': sqm_key,
    }
    if charge is None:
        cmd = '''export OMP_NUM_THREADS={thread};acpype -i {ligandfile} -b {molname} -a {forcefield} -c {method} -k "{sqm_key}" -f >acpype.log 2>&1 '''.format(**paras)
    else:
        cmd = '''export OMP_NUM_THREADS={thread};acpype -i {ligandfile} -b {molname} -a {forcefield} -c {method} -n {net_charge} -k "{sqm_key}" -f >acpype.log 2>&1 '''.format(**paras)
    RC = os.system(cmd)
    if RC != 0:
        print(cmd)
        os.system('tail -n 50 acpype.log')
        raise Exception('ERROR run the acpype. see the %s for details.'%os.path.abspath("acpype.log"))
    os.chdir(f'{molname}.acpype')
    moltop = pmd.load_file(f'{molname}_GMX.top')
    molitp = f'{molname}_GMX.itp'
    molgro = pmd.load_file(f'{molname}_GMX.gro', structure=True)
    if outtop:
        moltop.write(outtop)
    if outcoord:
        molgro.write_pdb(outcoord)
    if itpfile:
        RC = os.system(f'cp {molitp} {itpfile}')
    os.chdir(cwd)
    if not verbose:
        shutil.rmtree(ligandName)
    return moltop, molgro


def pdb2amber(pdbfile, outfile=None):
    structure = pmd.load_structure(pdbfile)
    atomlist = []
    for c, residue in enumerate(structure.residues, start=1):
        # change atoms name from GROMACS to AMBER
        for atom in residue.atoms:
            if atom.name == 'OC1':
                atom.name = 'O'
            elif atom.name == 'OC2':
                atom.name = 'OXT'
                residue.ter = True  # parmed terminal
        # change residues name according to AMBER
        if residue.name == 'ILE':
            for atom in residue.atoms:
                if atom.name == 'CD':
                    atom.name = 'CD1'
                    break
        elif residue.name == 'LYS':
            atoms = [atom.name for atom in residue.atoms]
            if 'HZ3' not in atoms:
                residue.name = 'LYN'
        elif residue.name == 'ASP':
            atoms = [atom.name for atom in residue.atoms]
            if 'HD2' in atoms:
                residue.name = 'ASH'
        elif residue.name == 'GLU':
            atoms = [atom.name for atom in residue.atoms]
            if 'HE2' in atoms:
                residue.name = 'GLH'
        elif residue.name in ['HIS', 'HIE', 'HID', 'HIP']:
            atoms = [atom.name for atom in residue.atoms if atom.atomic_number == 1]
            if 'HD1' in atoms and 'HE2' in atoms:
                residue.name = 'HIP'
            elif 'HD1' in atoms:
                residue.name = 'HID'
            elif 'HE2' in atoms:
                residue.name = 'HIE'
        for atom in residue.atoms:
            if atom.element > 2:
                atomlist.append(atom)
    structure.atoms = pmd.structure.AtomList(atomlist)
    if outfile is None:
        outfile = pdbfile
    structure.write_pdb(outfile)


def build_protein_tleap(pdbfile, forcefield='', outtop=None, outcoor=None):
    struc = pmd.load_structure(pdbfile)
    leaprc = '''
    source leaprc.protein.ff19SB
    source leaprc.gaff
    loadOff atomic_ions.lib
    loadamberparams frcmod.ions234lm_126_tip3p
    set default PBRadii mbondi2
    REC = loadpdb %s
    REC_OUT = combine { REC }
    saveamberparm REC_OUT REC.prmtop REC.inpcrd
    ''' % pdbfile


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
    uid = str(uuid.uuid4())
    proteinName = uid + os.path.split(pdbfile)[-1][:-4]+'.TOP'
    if not os.path.exists(proteinName):
        os.mkdir(proteinName)
    pdbfile = os.path.abspath(pdbfile)
    cwd = os.getcwd()
    os.chdir(proteinName)
    fix_insertions(pdbfile, 'begain.pdb')
    paras = {
        'gmx':GMXEXE,
        'pdbfile':'begain.pdb',
        'outfile': '1-pdb2gmx.gro',
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
    boxpdb = engine.gmx_box('1-pdb2gmx.gro', boxtype='triclinic', boxsize=0.9)
    #solpdb = engine.gmx_solvate(boxpdb, 'topol.top', maxsol=5)
    #ionspdb = engine.gmx_ions(solpdb, 'topol.top', conc=None, nNA=1, nCL=1, neutral=False)
    #engine.
    protgro = pmd.load_file('1-pdb2gmx.gro', structure=True)
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

def build_topol(receptor, ligand, outpdb, outtop, proteinforce='amber99sb-ildn', ligandforce='gaff2', charge_method='bcc', nt=1, verbose=False):
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
    indexfile = 'index.ndx'
    if isinstance(ligand, str):
        if ligand.endswith('.pdb'):
            moltop, molgro = build_protein(ligand, forcefield=proteinforce)
            gen_index_for_gbsa(prottop, moltop, indexfile)
        else:
            moltop, molgro = build_lignad(ligand, forcefield=ligandforce, charge_method=charge_method, nt=nt, verbose=verbose)
    elif ligand:
        moltop, molgro = ligand

    if ligand is None:
        systop, sysgro = prottop, protgro
    else:
        systop = moltop + prottop
        sysgro = molgro + protgro

    systop.write(outtop)
    sysgro.write_pdb(outpdb)
    newlines = []
    mtypes = []
    with open(outtop) as fr:
        lines = fr.readlines()
        for i,line in enumerate(lines):
            if line.startswith('[') and line.split()[1]=='system':
                fr = open(TEMPLATE_TOP)
                template = json.load(fr)
                fr.close()
                for k, v in template.items():
                    if k not in mtypes:
                        newlines.append(v)
            elif line.startswith('[') and line.split()[1]=='moleculetype':
                for l in lines[i+1:]:
                    if l.startswith('['):
                        break
                    elif l.strip() and not l.startswith(';'):
                        mtypes.append(l.split()[0])
                        break

            newlines.append(line)
    atomtypes = {
        'NA':  'Na      Na        22.990000  0.00000000  A        0.33284      0.0115897',
        'CL':  'Cl      Cl        35.450000  0.00000000  A       0.440104         0.4184',
        'SOL': 'OW      OW        16.000000  0.00000000  A       0.315061       0.636386\nHW      HW         1.008000  0.00000000  A              0              0',
    }
    with open(outtop, 'w') as fw:
        for line in newlines:
            fw.write(line)
            if line.startswith('[') and line.split()[1]=='atomtypes':
                for k,v in atomtypes.items():
                    if k not in mtypes:
                        fw.write(v+'\n')

    write_position_restrain(outtop)
    return indexfile

def main():
    pdbfile, ligandfile = sys.argv[1], sys.argv[2]

    grofile = 'sys.pdb'
    topfile = 'sys.top'

    build_topol(pdbfile, ligandfile, outpdb=grofile, outtop=topfile, nt=4)

    engine = GMXEngine()
    engine.run_to_md(grofile, topfile)

if __name__  == "__main__":
    main()
