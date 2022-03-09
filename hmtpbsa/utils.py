
import os
import shutil
from tokenize import group
def obtain_id_from_index(indexFile):
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

def generate_index_file(complexfile):
    
    cmd = '''gmx make_ndx -f %s 2>&1 << EOF
           name 2 LIGAND
           q
        EOF ''' % complexfile
    fr = os.popen(cmd)
    text = fr.read().strip()
    if 'Error' in text:
        print(cmd)
        raise Exception('Error run make_ndx: \n%s'%text)
    groupdict = {}
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
    groupdict['NACL'] = NACL
    cmd = '''gmx make_ndx -f {complexfile} -n index.ndx 2>&1 << EOF
        {non-Water} & {NACL} & ! {LIGAND}
        name {RECEPTOR} RECEPTOR
           q
        EOF '''.format(**groupdict)
    fr = os.popen(cmd)
    text = fr.read().strip()
    if 'Error' in text:
        print(cmd)
        raise Exception('Error run make_ndx: \n%s'%text)
    #shutil.rmtree('#index.ndx.1#')
    indexfile = os.path.abspath('index.ndx')
    return indexfile