import os
import pandas as pd
from turtle import title
from typing import List
from hmtpbsa.settings import GMXEXE
def set_amber_home(proc):
    """
    Find the directory containing the executable for a command
  
    Args:
      proc: The name of the executable to find.
  
    Returns:
      the path to the amber home directory.
    """
    cmd = 'which %s '%proc
    f = os.popen(cmd)
    text = f.read().strip()
    if not text:
        raise Exception("Command not found: %s "%proc)
    bindir = os.path.split(text)[0]
    amberhome = os.path.split(bindir)[0]
    return amberhome

def obtain_num_of_frame(trajfile):
    """
    Get the number of frames in a trajectory file
    
    Args:
      trajfile: the trajectory file, in .xtc or .trr format.
    
    Returns:
      the number of frames in the trajectory file.
    """
    cmd = '%s check -f %s 2>&1 |grep Coords'% (GMXEXE, trajfile)
    fr = os.popen(cmd)
    text = fr.read().strip()
    if not text:
        print(cmd)
        raise Exception("ERROR obtain %s's frame number.")
    nframe = int(text.split()[1])
    return nframe

def read_MMPBSA_output(datfile):
    DeltaG = {}
    datalist = []
    with open(datfile) as fr:
        for line in fr:
            line = line.strip()
            if line.startswith('|') or line.startswith('-') or not line:
                continue
            if line.startswith('GENERALIZED BORN'):
                tagName = 'GB'
            elif line.startswith('POISSON BOLTZMANN'):
                tagName = 'PB'
            elif ':' in line:
                groupname = line[:-1]
            else:
                component = line[:15].strip()
                Llist = line[15:].strip().split()
                if len(Llist)==5:
                    # method, group name, Energy Component, Average, SD(Prop.), SD, SEM(Prop.) SEM
                    datalist.append( (tagName, groupname, component, float(Llist[0]), float(Llist[1]), float(Llist[2]), float(Llist[3]), float(Llist[4])))
    title = ['type', 'group', 'component', 'average', 'SD(Prop.)', 'SD', 'SEM(Prop.)', 'SEM']
    df = pd.DataFrame(datalist, columns=title)
    return df

def read_DECOMP_output(datfile):
    dic = {}   
    model = {
        'Generalized Born model':'GB',
        'Poisson Boltzmann model':'PB'
    }
    header = None
    data = []
    with open(datfile) as fr:
        lines = fr.readlines()
    for i,line in enumerate(lines):
        line = line.strip()
        if line.startswith('Energy Decomposition Analysis'):
            if etype is not None and len(data)>0:
                dic[etype] = pd.DataFrame(data, columns=header)
                data = []
            etype = line.split(":")[1].strip()
            etype = model[etype]
            if header is None:
                tmp = line[i+4].split(',')
                header = ['residue']
                for h in tmp[1:]:
                    if h.strip():
                        header.append([h, h+' Std', h+' Std. Err. of Mean'])
        elif line.startswith(('L:', 'R:')):
            tmp = tmp.split(',')
            data.append(tmp)
    if len(data)>0:
        dic[etype] = pd.DataFrame(data, columns=header)
    return dic


        
