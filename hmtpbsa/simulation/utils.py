
import os

def convert_to_mol2(inputfile, filetype, outfile=None):
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
        outfile = filename + '.mol2'
    # convert to mol2
    cmd = 'obabel -i %s %s -o mol2 -O %s'%(filetype, inputfile, outfile)
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
    cmd = '''gmx make_ndx -f %s 2>&1 << EOF
           name 2 LIGAND
           q
        EOF ''' % complexfile
    fr = os.popen(cmd)
    text = fr.read().strip()
    if 'Error' in text:
        print(cmd)
        raise Exception('Error run make_ndx: \n%s'%text)
    groupdict = {
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
    cmd = '''gmx make_ndx -f {complexfile} -n index.ndx 2>&1 <<EOF
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