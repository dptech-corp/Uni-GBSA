import os
import uuid
import shutil
from pathlib import Path
from unigbsa.settings import GMXEXE, logging


def convert_format(inputfile, filetype, outfile=None, outtype='mol'):
    """
    Convert a file of type filetype to mol2 format
    
    Args:
      inputfile: the input file name
      filetype: the type of input file.
      outfile: the name of the output file. If not specified, the output file name will be the same as
    the input file name, but with a different extension.
    
    Returns:
      the name of the output file.
    """
    if outfile is None:
        filename = os.path.split(inputfile)[-1][:-4]
        outfile = filename + '.' + outtype
    # convert to mol2
    cmd = 'obabel -i %s %s -o %s -O %s >/dev/null 2>&1 '%(filetype, inputfile, outtype, outfile)
    RC = os.system(cmd)
    if RC!=0:
        raise Exception('ERROR: failed convert %s to %s'%(inputfile, outfile))
    return os.path.abspath(outfile)

def assign_partial_charge(inputfile, filetype, charge_method='gasteiger', outfile=None):
    '''

eem    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Bultinck B3LYP/6-31G*/MPA
eem2015ba    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf B3LYP/6-311G/AIM
eem2015bm    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf B3LYP/6-311G/MPA
eem2015bn    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf B3LYP/6-311G/NPA
eem2015ha    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf HF/6-311G/AIM
eem2015hm    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf HF/6-311G/MPA
eem2015hn    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf HF/6-311G/NPA
eqeq    Assign EQEq (charge equilibration) partial charges.
fromfile    Assign charges from file containing {'atom-name', charge} pairs
gasteiger    Assign Gasteiger-Marsili sigma partial charges
mmff94       Assign MMFF94 partial charges
none    Clear all partial charges
qeq    Assign QEq (charge equilibration) partial charges (Rappe and Goddard, 1991)
qtpie    Assign QTPIE (charge transfer, polarization and equilibration) partial charges (Chen and Martinez, 2007)
    '''
    charge_methods = ['eem', 'eem2015ba', 'eem2015bm', 'eem2015bn', 'eem2015ha', 'eem2015hm', 'eem2015hn', 'gasteiger', 'mmff94', 'qeq', 'qtpie']
    if outfile is None:
        filename = os.path.split(inputfile)[-1][:-4]
        outfile = filename + '.mol2'
    if charge_method not in charge_methods:
        raise Exception('ERROR: charge method %s is not one of the %s '%(charge_method, ','.join(charge_methods)))
    # convert to mol2
    cmd = 'obabel -i %s %s -o %s -O %s -xl --partialcharge %s >/dev/null 2>&1 '%(filetype, inputfile, 'mol2', outfile, charge_method)
    RC = os.system(cmd)
    if RC!=0:
        raise Exception('ERROR: failed convert %s to %s'%(inputfile, outfile))
    return os.path.abspath(outfile)

def guess_filetype(inputfile):
    """Guess the file type of input file.

    Args:
        inputfile ([str]): A input file to guess filetype

    Returns:
        [str]: A filetype
    """
    basefile = os.path.split(inputfile)[-1]
    filetype = basefile.split('.')[-1]
    return filetype

def load_position_restraints(topfile, outfile=None):
    """
    Reads a topology file and writes a new topology file with position restraints included
    
    Args:
      topfile: the topology file
      outfile: the output file name.
    """
    if outfile is None:
        outfile = topfile
    with open(topfile) as fr:
        lines = fr.readlines()
    POSRES = False
    fw = open(outfile, 'w')
    for line in lines:
        if line.strip().endswith('#ifdef POSRES'):
            POSRES = True
        if line.startswith('#include') and POSRES:
            itpfile = line.split()[1].strip().replace('"','')
            fw.write("\n")
            print(itpfile)
            with open(itpfile) as fr:
                for line in fr:
                    if not line.startswith(';'):
                        fw.write(line)
            fw.write("\n")
            POSRES = False
        else:
            fw.write(line)

def generate_index_file_for_restrain(complexfile):
    cmd = '''%s make_ndx -f %s 2>&1 << EOF
           name 2 LIGAND
           q
        EOF ''' % ( GMXEXE, complexfile )
    fr = os.popen(cmd)
    text = fr.read().strip()
    if 'Error' in text:
        print(cmd)
        raise Exception('Error run make_ndx: \n%s'%text)
    groupdict = {
        'gmx':GMXEXE
        }
    groupid = 0
    with open('index.ndx') as fr:
        for line in fr:
            if line.startswith('['):
                tmp = line.split()
                groupdict[tmp[1]] = groupid
                groupid += 1
    groupdict['RECEPTOR'] = groupid
    groupdict['H_atoms'] = groupid + 1
    groupdict['all_heavy'] = groupid + 2
    groupdict['complex_heavy'] = groupid + 3
    groupdict['complexfile'] = complexfile
    NACL = ''
    if 'NA' in groupdict:
        NACL += ' ! %d &'%groupdict['NA']
    if 'CL' in groupdict:
        NACL += ' ! %d'%groupdict['CL']
    if 'non-Water' not in groupdict:
        groupdict['non-Water'] = ''

    groupdict['NACL'] = NACL
    cmd = '''{gmx} make_ndx -f {complexfile} -n index.ndx 2>&1 <<EOF
        ! {LIGAND} & {NACL} & {non-Water}
        name {RECEPTOR} RECEPTOR
        {NACL} & a H*
        name {H_atoms} H_atoms
        ! a H*
        name {all_heavy} all_heavy
        {LIGAND} | {RECEPTOR} & ! {H_atoms}
        name {complex_heavy} complex_heavy

           q\nEOF'''.format(**groupdict)
    #print(cmd)
    fr = os.popen(cmd)
    text = fr.read().strip()
    if 'Error' in text:
        print(cmd)
        raise Exception('Error run make_ndx: \n%s'%text)
    #shutil.rmtree('#index.ndx.1#')
    os.system('rm \#*')
    indexfile = os.path.abspath('index.ndx')
    return indexfile

def write_position_restrain(topfile, outfile=None, fc=[1000, 1000, 1000], excludes=['NA','CL']):
    if outfile is None:
        outfile = topfile
    with open(topfile) as fr:
        lines = fr.readlines()
    watersName = ['SOL', 'TIP3P']
    fw = open(outfile, 'w')
    positionRestrain = []
    moleculeName = None
    linetype = None
    for i, line in enumerate(lines):
        if line.strip().startswith('['):
            linetype = line.split('[')[-1].split(']')[0].strip()
            if linetype == 'moleculetype' and len(positionRestrain)>0:
                if moleculeName in watersName:
                    fw.write('#ifdef POSRES_H2O\n[ position_restraints ]\n;  i funct       fcx        fcy        fcz\n')
                else:
                    fw.write('#ifdef POSRES_HEAVY\n[ position_restraints ]\n;  i funct       fcx        fcy        fcz\n')
                fw.writelines(positionRestrain)
                fw.write('#endif\n\n')
                positionRestrain = []
        elif line.strip().startswith('#'):
            linetype = 'define'
        elif line.strip() and not line.strip().startswith(';'):
            if linetype == 'moleculetype':
                    moleculeName = line.strip().split()[0]
                    if moleculeName in excludes:
                        moleculeName = None
            elif moleculeName and linetype == 'atoms':
                lineList = line.strip().split()
                atomindex, atomtype = lineList[0], lineList[1].upper()
                if not atomtype.startswith('H'):
                    positionRestrain.append('%5s   1   %5d %5d %5d\n'%(atomindex, fc[0], fc[1], fc[2]))
        fw.write(line)
    fw.close()
    return outfile

def fix_insertions(pdbfile, outfile=None):
    """
    Delete insertion codes (at specific residues).

    By default, removes ALL insertion codes on ALL residues. Also bumps
    the residue numbering of residues downstream of each insertion.

    This function is a generator.

    Parameters
    ----------
    fhandle : a line-by-line iterator of the original PDB file.

    option_list : list
        List of insertion options to act on.
        Example ["A9", "B12"]. An empty list performs the default
        action.

    Yields
    ------
    str (line-by-line)
        The modified (or not) PDB line.
    """

    # Keep track of residue numbering
    # Keep track of residues read (chain, resname, resid)
    offset = 0
    prev_resi = None
    seen_ids = set()
    clean_icode = False
    records = ('ATOM', 'HETATM', 'ANISOU')
    with open(pdbfile) as fr:
        lines = fr.readlines()
    if outfile is None:
        outfile = pdbfile
    fw = open(outfile, 'w')
    residMapping = {}
    residMappingBack = {}
    has_insertion = False
    for line in lines:

        if line.startswith(records):
            res_uid = line[17:27]  # resname, chain, resid, icode
            id_res = line[21] + line[22:26].strip()  # A99, B12
            has_icode = line[26].strip()  # ignore ' ' here
            # unfortunately, this is messy but not all PDB files follow a nice
            # order of ' ', 'A', 'B', ... when it comes to insertion codes..
            if prev_resi != res_uid:  # new residue
                # Does it have an insertion code
                # OR have we seen this chain + resid combination before?
                # #2 to catch insertions WITHOUT icode ('A' ... ' ' ... 'B')
                if (has_icode or id_res in seen_ids):
                    # Do something about it
                    # if the user provided options and this residue is in them
                    # OR if the user did not provide options
                    clean_icode = True
                else:
                    clean_icode = False

                prev_resi = res_uid

                if id_res in seen_ids:  # offset only if we have seen this res.
                    offset += 1

            # Modify resid if necessary
            oldresid = int(line[22:26])
            resid = oldresid + offset
            chainID = line[21]
            if chainID not in residMapping:
                residMapping[chainID] = {str(oldresid):resid}
                residMappingBack[chainID] = {resid:str(oldresid)}
            else:
                if str(oldresid) in residMapping[chainID]:
                    residMapping[chainID].update({"%d%s"%(oldresid, line[26].strip()):resid})
                    residMappingBack[chainID].update({resid:"%d%s"%(oldresid, line[26].strip())})
                else:
                    residMapping[chainID].update({str(oldresid):resid})
                    residMappingBack[chainID].update({resid:str(oldresid)})
            if clean_icode:  # remove icode
                line = line[:26] + ' ' + line[27:]
            line = line[:22] + str(resid).rjust(4) + line[26:]
            seen_ids.add(id_res)

            # Reset offset on TER
            if line.startswith('TER'):
                offset = 0

        fw.write(line)
    fw.close()
    return residMapping, residMappingBack


def prepare_ligand(molfile, outfile=None):
    from rdkit import Chem
    if outfile is None:
        outfile = molfile
    mol = Chem.MolFromMolFile(molfile)
    mol_out = Chem.rdmolops.AddHs(mol, addCoords=True)
    Chem.MolToMolFile(mol_out, outfile)
    return outfile


def obtain_net_charge_command(molfile):
    import uuid
    ftype = guess_filetype(molfile)
    mol2file = str(uuid.uuid4()) + '.mol2'
    cmd = f'obabel -i {ftype} {molfile} -omol2 -O {mol2file}  >/dev/null 2>&1 '
    RC = os.system(cmd)
    if RC != 0:
        print('Failed to obtain mol charge, use guess')
        return None
    charge = 0
    with open(mol2file) as fr:
        for line in fr:
            llist = line.split()
            if len(llist) >= 8:
                charge += float(llist[-1])
    return int(round(charge))


def obtain_net_charge_rdkit(sdfile):
    from rdkit import Chem
    # Read the MOL or SDF file
    mol = Chem.MolFromMolFile(sdfile)

    # Calculate the net charge of the molecule
    net_charge = Chem.GetFormalCharge(mol)

    return net_charge


def obtain_net_charge(sdfile):
    from openbabel import openbabel
    # Read the MOL or SDF file
    mol = openbabel.OBMol()
    with open(sdfile, 'r') as file:
        format = sdfile[-3:]
        file.seek(0)  # Reset the file pointer
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat(format)
        obConversion.ReadString(mol, file.read())

    # Calculate the net charge of the molecule
    net_charge = mol.GetTotalCharge()
    return net_charge


def check_forcefield(sdfile):
    sqm_key = "grms_tol=0.005,qm_theory='AM1',scfconv=1.d-5,ndiis_attempts=700,maxcyc=0"
    tmpdir = Path('/tmp') / str(uuid.uuid4())
    sdfile = Path(sdfile).absolute()
    net_charge = obtain_net_charge(str(sdfile))
    tmpdir.mkdir(exist_ok=True)
    cwd = os.getcwd()
    os.chdir(tmpdir)
    cmd = f'export OMP_NUM_THREADS=1;acpype -i {sdfile} -b MOL -a gaff -c bcc -n {net_charge} -k "{sqm_key}" >acpype.log 2>&1 '
    rc = os.system(cmd)
    os.chdir(cwd)
    shutil.rmtree(tmpdir)
    if rc != 0:
        return False
    else:
        return True


def check_element(sdfile):
    from openbabel import openbabel
    ext = Path(sdfile).suffix[1:]
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(ext, ext)
    mol = openbabel.OBMol()
    if not obConversion.ReadFile(mol, str(sdfile)):
        return -1

    acceptable_elements = {6, 7, 8, 16, 15, 1, 9, 17, 35, 53}
    for atom in openbabel.OBMolAtomIter(mol):
        element = atom.GetAtomicNum()
        if element not in acceptable_elements:
            return 1
    return 0


def get_total_valence_electrons(sdfile):
    from openbabel import openbabel
    extin = Path(sdfile).suffix[1:]
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(extin, extin)
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, str(sdfile))
    total_valence_electrons = 0
    for atom in openbabel.OBMolAtomIter(mol):
        total_valence_electrons += atom.GetAtomicNum() - atom.GetFormalCharge()
    return total_valence_electrons


def get_electronegativity(atomic_number):
    electronegativity_table = {
        1: 2.20,  # Hydrogen (H)
        6: 2.55,  # Carbon (C)
        7: 3.04,  # Nitrogen (N)
        8: 3.44,  # Oxygen (O)
        15: 2.19,  # Phosphorus (P)
        16: 2.58,  # Sulfur (S)
        9: 3.98,  # Fluorine (F)
        17: 3.16,  # Chlorine (Cl)
        35: 2.96,  # Bromine (Br)
        53: 2.66,  # Iodine (I)
    }
    return electronegativity_table.get(atomic_number, 0)


def adjust_charge_based_on_electronegativity(ob_mol):
    from openbabel import openbabel
    max_electronegativity = 0
    target_atom = None

    for atom in openbabel.OBMolAtomIter(ob_mol):
        if atom.GetFormalCharge() != 0:
            electronegativity = get_electronegativity(atom.GetAtomicNum())
            if electronegativity > max_electronegativity:
                max_electronegativity = electronegativity
                target_atom = atom

    if target_atom is not None:
        if target_atom.GetFormalCharge() > 0:
            target_atom.SetFormalCharge(target_atom.GetFormalCharge() - 1)
        elif target_atom.GetFormalCharge() < 0:
            target_atom.SetFormalCharge(target_atom.GetFormalCharge() + 1)


def add_hydrogen(sdfile, outfile=None):
    from openbabel import openbabel
    extin = extout = Path(sdfile).suffix[1:]
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(extin, extout)
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, str(sdfile))
    mol.AddHydrogens()
    mol_string = obConversion.WriteString(mol)
    if outfile:
        with open(outfile, 'w') as fw:
            fw.write(mol_string)
    else:
        return mol


def repair_ligand(sdfile, outfile=None):
    from openbabel import openbabel
    extin = Path(sdfile).suffix[1:]
    if outfile is None:
        fname = str(uuid.uuid4()) + '.' + extin
        outfile = Path('/tmp') / fname
    extout = Path(outfile).suffix[1:]
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(extin, extout)
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, str(sdfile))
    mol.AddHydrogens()
    mol.DeleteHydrogens()
    # Correct valence information
    pH = 7.4
    charge_model = openbabel.OBChargeModel.FindType("mmff94")
    charge_model.ComputeCharges(mol, str(pH))
    mol.CorrectForPH()
    mol.PerceiveBondOrders()
    total_valence_electrons = sum([atom.GetAtomicNum() - atom.GetFormalCharge() for atom in openbabel.OBMolAtomIter(mol)])
    if total_valence_electrons % 2 == 1:
        adjust_charge_based_on_electronegativity(mol)
    mol_string = obConversion.WriteString(mol)
    new_string = []
    for i, line in enumerate(mol_string.split('\n')):
        if i >= 4 and len(line) >= 30 and line.startswith(' '):
            line = line[:42] + '  0  0  0  0  0  0  0  0  0'
        new_string.append(line)
    mol_string = '\n'.join(new_string)
    with open(outfile, 'w') as fw:
        fw.write(mol_string)
    add_hydrogen(outfile, outfile)
    return outfile


def ligand_validate(sdfile, outfile=None):
    rc = check_element(sdfile)
    if rc == -1:
        logging.error(f'Failed to load {sdfile}, please check your input {sdfile}.')
        exit(256)
    elif rc == 1:
        logging.error(f'Ligand file only accept C N O S P H F Cl Br I, please check your input {sdfile}.')
        exit(257)
    total_valence_electrons = get_total_valence_electrons(sdfile)
    if total_valence_electrons % 2 == 1 or not check_forcefield(sdfile):
        logging.warning(f'The total valence electrons of your ligand is odd({total_valence_electrons}) or forcefield check error, try to repair input ligand.')
        outfile = repair_ligand(sdfile, outfile=outfile)
        total_valence_electrons = get_total_valence_electrons(outfile)
        if total_valence_electrons % 2 == 1 or not check_forcefield(outfile):
            logging.error(f'The total valence electrons of your ligand is stil odd({total_valence_electrons}) or forcefield check error after repair ligand, please check your input {sdfile}.')
            exit(257)
        else:
            return outfile
    else:
        return sdfile


def gen_index_for_gbsa(rec, lig, outfile='index.ndx'):
    lig_atoms, rec_atoms = 0, 0
    for key, (mol, _) in rec.molecules.items():
        if key.startswith(('MOL', 'protein', 'system', 'Protein')):
            rec_atoms += len(mol.atoms)
    for key, (mol, _) in lig.molecules.items():
        if key.startswith(('MOL', 'protein', 'system', 'Protein')):
            lig_atoms += len(mol.atoms)

    atom_dict = {
        'System': (1, lig_atoms+rec_atoms+1),
        'ligand': (1, lig_atoms+1),
        'receptor': (lig_atoms+1, rec_atoms+lig_atoms+1)
    }
    with open(outfile, 'w') as fw:
        for k, (s, e) in atom_dict.items():
            fw.write(f'\n[ {k} ]\n')
            for i, ndx in enumerate(range(s, e)):
                fw.write(str(ndx) + ' ')
                if i != 0 and i % 10 == 0:
                    fw.write('\n')
