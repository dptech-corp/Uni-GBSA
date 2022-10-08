from email.policy import default
import os
import shutil
import configparser
from unigbsa.settings import logging, GMXEXE, DEFAULT_CONFIGURE_FILE

def obtain_id_from_index(indexFile):
    """
    Given an index file, return the group IDs of the receptor and ligand
    
    Args:
      indexFile: the index file that contains the information about the receptor and ligand groups.
    
    Returns:
      the receptor and ligand IDs.
    """
    receptorID = None
    ligandID = None
    groupCount = 0
    groupNames = []
    with open(indexFile) as fr:
        for line in fr:
            if line.startswith('['):
                groupName = line.strip()[1:-1].strip()
                if groupName.lower() == 'receptor':
                    receptorID = groupCount
                elif groupName.lower() == 'ligand':
                    ligandID = groupCount
                groupCount += 1
                groupNames.append(groupName)
    if receptorID is None:
        print(groupNames)
        raise Exception('Not found "receptor" group in the index file: %s'%indexFile)
    if ligandID is None:
        print(groupNames)
        raise Exception('Not found "ligand" group in the index file: %s'%indexFile)
    return receptorID, ligandID

def generate_index_file(complexfile, pbc=False):
    """
    Generate index file for the complex file
    
    Args:
      complexfile: The name of the complex file.
    
    Returns:
      The index file is being returned.
    """
    
    cmd = '''%s make_ndx -f %s 2>&1 << EOF
           name 2 LIGAND
           q
        EOF''' % (GMXEXE, complexfile)
    fr = os.popen(cmd)
    text = fr.read().strip()
    if 'Error' in text:
        print(cmd)
        raise Exception('Error run make_ndx: \n%s'%text)
    groupdict = {
        'center': '',
        'output': '',
        'gmx':GMXEXE,
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
    groupdict['complexfile'] = complexfile
    NACL = ''
    if 'NA' in groupdict:
        NACL += ' ! %d &'%groupdict['NA']
    if 'CL' in groupdict:
        NACL += ' ! %d'%groupdict['CL']
    if 'non-Water' not in groupdict:
        groupdict['non-Water'] = ''
    if pbc:
        groupdict['center'] = '2 \n name %d center'%(groupid+1)
        groupdict['output'] = '2|%d \n name %d output'%(groupdict['RECEPTOR'], groupid+2)
    groupdict['NACL'] = NACL
    cmd = '''{gmx} make_ndx -f {complexfile} -n index.ndx 2>&1 <<EOF
        ! {LIGAND} & {NACL} & {non-Water}
        name {RECEPTOR} RECEPTOR
        {center}
        {output}
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

def process_pbc(trajfile, tprfile, indexfile, outfile=None, logfile="/dev/null"):
    fname = os.path.split(trajfile)[-1][:-4]
    suffix = os.path.split(trajfile)[-1][-4:]
    if outfile is None:
        outfile = '%s-pbc'%fname + suffix
    basecmd = '%s trjconv -s %s -n %s '%(GMXEXE, tprfile, indexfile)
    cleanfiles = []
    ## step 1:
    inf  = trajfile
    outf = '%s-1'%fname + suffix
    cmd = basecmd + '''-f %s -o %s -pbc mol -ur compact >> %s 2>&1 << EOF
    output\nEOF'''%(inf, outf, logfile)
    RC = os.system(cmd)
    if RC != 0:
        print(cmd)
        raise Exception('Process pbc error: see the %s file for details'%logfile)
    cleanfiles.append(outf)
    ## step 2: center
    inf = outf
    outf = '%s-2'%fname + suffix
    cmd = basecmd + '''-f %s -o %s -center >> %s 2>&1 <<EOF
    center
    output\nEOF'''%(inf, outf, logfile)
    RC = os.system(cmd)
    if RC != 0:
        print(cmd)
        raise Exception('Process pbc error: see the %s file for details'%logfile)
    cleanfiles.append(outf)
    ## step3:
    inf = outf
    outf = '%s-3'%fname + suffix
    cmd = basecmd + '''-f %s -o %s -pbc res -ur compact >> %s 2>&1 <<EOF
    output\nEOF'''%(inf, outf, logfile)
    RC = os.system(cmd)
    if RC != 0:
        print(cmd)
        raise Exception('Process pbc error: see the %s file for details'%logfile)
    cleanfiles.append(outf)
    ## step4:
    inf = outf
    outf = '%s-4'%fname + suffix
    cmd = basecmd + '''-f %s -o %s -pbc mol -ur compact >> %s 2>&1 <<EOF
    output\nEOF'''%(inf, outf, logfile)
    RC = os.system(cmd)
    if RC != 0:
        print(cmd)
        raise Exception('Process pbc error: see the %s file for details'%logfile)
    os.system('mv %s %s'%(outf, outfile))

    ## step5: remove
    cleanline = ' '.join(cleanfiles)
    cmd = "rm -rf %s"%cleanline
    os.system(cmd)

    return outfile

def load_configue_file(conf=None):
    import json
    if conf is None:
        conf = DEFAULT_CONFIGURE_FILE
    fileformat = os.path.split(conf)[-1].split('.')[-1]
    if fileformat == 'json':
        with open(conf) as fr:
            paras =  json.load(fr)
        default = {
                "mode": "em",
                "boxtype": "triclinic",
                "boxsize": 0.9,
                "conc": 0.15,
                "nsteps": 500000,
                "nframe": 100,
                "eqsteps": 50000,
                "proteinforcefield": "amber03",
                "ligandforcefield": "gaff",
                "maxsol": 0,
                "ligandCharge": "bcc"
            }
        for k,v in default.items():
            if k not in paras['simulation']:
                paras['simulation'][k] = v
    elif fileformat == 'ini':
        config = configparser.ConfigParser()
        config.read(conf)
        paras = {
                'simulation':{
                'mode': config.get('simulation', 'mode', fallback='em'),
                'boxtype' : config.get('simulation', 'boxtype', fallback='triclinic'),
                'boxsize' : config.getfloat('simulation', 'boxsize', fallback=0.9),
                'conc': config.getfloat('simulation', 'conc', fallback=0.15),
                'nsteps': config.getint('simulation', 'nsteps', fallback=500000),
                'nframe': config.getint('simulation', 'nframe', fallback=100),
                'eqsteps': config.getint('simulation', 'eqsteps', fallback=50000),
                'proteinforcefield': config.get('simulation', 'proteinforcefield', fallback='amber99sb-ildn'),
                'ligandforcefield': config.get('simulation', 'ligandforcefield', fallback='gaff2'),
                'maxsol': config.getint('simulation', 'maxsol', fallback=0),
                'ligandCharge': config.get('simulation', 'ligandCharge', fallback='bcc'),
            },
            'GBSA':  {k:v for k,v in config.items('GBSA')}
        }

    return paras