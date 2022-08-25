
import pandas as pd
def read_FINAL_output(datfile):
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
            if len(data)>0:
                dic[etype] = pd.DataFrame(data, columns=header)
                data = []
            tagName = None
            etype = line.split(":")[1].strip()
            etype = model[etype]
            if header is None:
                tmp = lines[i+4].split(',')
                header = ['residue']
                for h in tmp[1:]:
                    if h.strip():
                        header.extend([h, h+' Std', h+' Std. Err. of Mean'])
        elif line.startswith(('L:', 'R:')) and tagName=='DELTAS:':
            tmp = line.split(',')
            data.append(tmp)
        elif line.startswith('Total Energy Decomposition'):
            tagName = lines[i-1].strip()
        elif line.startswith('Residue') and not lines[i-1].startswith('Total Energy Decomposition'):
            tagName = ''
    if len(data)>0:
        dic[etype] = pd.DataFrame(data, columns=header)
    return dic

def read_DEO_output(csvfile):
    dic = {}
    header = None
    data = []
    dtype = None
    with open(csvfile) as fr:
        lines = fr.readlines()
    for i,line in enumerate(lines):
        if "Poisson Boltzmann" in line:
            if len(data)!=0:
                dic[ptype] = pd.DataFrame(data, columns=header)
                data = []
            ptype = 'PB'
        elif "Generalized Born" in line:
            if len(data)!=0:
                dic[ptype] = pd.DataFrame(data, columns=header)
            ptype = 'GB'
        elif line.startswith('Total Decomposition Contribution'):
            dtype = lines[i-1].strip()
            if header is None:
                header = lines[i+1].strip().split(',')
        elif dtype == 'DELTAS:' and ('L:' in line or 'R:' in line):
            data.append(line.strip().split(','))
    if len(data)!=0:
        dic[ptype] = pd.DataFrame(data, columns=header)
    return dic

def read_EO_output(csvfile):
    dic = {}
    data = []
    header = None
    dtype = None
    ptype = None
    with open(csvfile) as fr:
        lines = fr.readlines()
    for i,line in enumerate(lines):
        if "POISSON BOLTZMANN" in line:
            if len(data)!=0:
                dic[ptype] = pd.DataFrame(data, columns=header)
                data = []
            ptype = 'PB'
        elif "GENERALIZED BORN" in line:
            if len(data)!=0:
                dic[ptype] = pd.DataFrame(data, columns=header)
                data = []
            ptype = 'GB'
        elif line.startswith('Frame #'):
            dtype = lines[i-1].strip()
            header = lines[i].strip().split(',')
        elif dtype == 'Delta Energy Terms' and line.startswith(tuple(map(str, range(10)))):
            data.append(line.strip().split(','))
    if len(data)!=0:
        dic[ptype] = pd.DataFrame(data, columns=header)
    return dic