

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