
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
    return outfile

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

