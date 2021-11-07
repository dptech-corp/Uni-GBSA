import argparse
import os
import sys
from xml.sax import parse
import numpy as np
import pandas as pd
import MDAnalysis as mda
import logging

ProcPath = os.path.abspath(os.path.dirname(sys.argv[0]))
__ResidueTypeFile = os.path.join(ProcPath, 'data/ResidueType.pkl3')
ResnameDictionary = pd.read_pickle(__ResidueTypeFile)


def detect_group(pdbfile:str) -> dict:
    def map_restype(resname:str):
        if resname in ResnameDictionary:
            ResType = ResnameDictionary[resname]
        else:
            ResType = 'Other'
        return ResType

    u = mda.Universe(pdbfile)
    if hasattr(u.residues, 'record_types'):
        recordTypes = [ recod[0] for recod in u.residues.record_types ]
    else:
        recordTypes = ['ATOM'] * u.residues.n_residues
    
    resids = u.residues.resids
    resnames = u.residues.resnames
    chainIDs = u.residues.segids
    Restype = list(map(map_restype, resnames))
    residsContinue = [0]
    index = 0

    for i in range(1, len(resids)):
        if ((recordTypes[i] == 'HETATM' or Restype[i] in ['ion','drug','solvent','ion complex','other', 'coenzyme'])\
              and resnames[i] != resnames[i-1]
            ) \
          or (resids[i] - resids[i-1]) != 1:
            index += 1
        residsContinue.append(index)
    dic = {
        'recordTypes':recordTypes,
        'resids':resids,
        'resnames':resnames,
        'chain':chainIDs,
        'restype': Restype,
        'continue':residsContinue
    }
    df = pd.DataFrame(dic)
    groupKeys = ['recordTypes', 'chain', 'continue']#, 'restype']
    groups = df.groupby(by=groupKeys)

    GroupList = []
    for _, group in groups:
        recordType = group['recordTypes'].values[0]
        ChainId = group['chain'].values[0]
        resnames = group['resnames'].values
        ResType = group['restype'].values[0]
        resid = group['resids'].values
        sortKey = '%s%d'%(ChainId, resid[0]+10000000)
        GroupList.append((index, recordType, ChainId, resid, resnames, ResType, sortKey))
    GroupDict = {}
    for index, group in enumerate(sorted(GroupList, key= lambda x: x[6])):
        RecordType = group[1]
        ChainId = group[2]
        ResidRange = str(group[3][0]) if len(group[3])==1 else "%d:%d"%(group[3][0], group[3][-1])
        resnames = group[4]
        ResType = group[5]
        GroupDict[index] = (index, RecordType, ChainId, ResidRange, resnames, ResType)

    return GroupDict


def detect_group_old(pdbfile: str) -> dict:
    """Detect the group information from the input pdb file.

    Args:
        pdbfile (str): input pdb format file.

    Returns:
        dict: A dict of group infomation.
        example:
        {
            0: (0, 'ATOM', 'A', '2:306', ['ASP', 'THR', 'ARG', 'ILE', 'GLY', 'VAL', 'THR', 'ILE', 'TYR', 'LYS', 'TYR', 'ASP', 'ASP', 'ASN', 'PHE', 'MET', 'SER', 'VAL', 'VAL', 'ARG', 'LYS', 'ALA', 'ILE', 'GLU', 'GLN', 'ASP', 'ALA', 'LYS', 'ALA', 'ALA', 'PRO', 'ASP', 'VAL', 'GLN', 'LEU', 'LEU', 'MET', 'ASN', 'ASP', 'SER', 'GLN', 'ASN', 'ASP', 'GLN', 'SER', 'LYS', 'GLN', 'ASN', 'ASP', 'GLN', 'ILE', 'ASP', 'VAL', 'LEU', 'LEU', 'ALA', 'LYS', 'GLY', 'VAL', 'LYS', 'ALA', 'LEU', 'ALA', 'ILE', 'ASN', 'LEU', 'VAL', 'ASP', 'PRO', 'ALA', 'ALA', 'ALA', 'GLY', 'THR', 'VAL', 'ILE', 'GLU', 'LYS', 'ALA', 'ARG', 'GLY', 'GLN', 'ASN', 'VAL', 'PRO', 'VAL', 'VAL', 'PHE', 'PHE', 'ASN', 'LYS', 'GLU', 'PRO', 'SER', 'ARG', 'LYS', 'ALA', 'LEU', 'ASP', 'SER', 'TYR', 'ASP', 'LYS', 'ALA', 'TYR', 'TYR', 'VAL', 'GLY', 'THR', 'ASP', 'SER', 'LYS', 'GLU', 'SER', 'GLY', 'ILE', 'ILE', 'GLN', 'GLY', 'ASP', 'LEU', 'ILE', 'ALA', 'LYS', 'HIS', 'TRP', 'ALA', 'ALA', 'ASN', 'GLN', 'GLY', 'TRP', 'ASP', 'LEU', 'ASN', 'LYS', 'ASP', 'GLY', 'GLN', 'ILE', 'GLN', 'PHE', 'VAL', 'LEU', 'LEU', 'LYS', 'GLY', 'GLU', 'PRO', 'GLY', 'HIS', 'PRO', 'ASP', 'ALA', 'GLU', 'ALA', 'ARG', 'THR', 'THR', 'TYR', 'VAL', 'ILE', 'LYS', 'GLU', 'LEU', 'ASN', 'ASP', 'LYS', 'GLY', 'ILE', 'LYS', 'THR', 'GLU', 'GLN', 'LEU', 'GLN', 'LEU', 'ASP', 'THR', 'ALA', 'MET', 'TRP', 'ASP', 'THR', 'ALA', 'GLN', 'ALA', 'LYS', 'ASP', 'LYS', 'MET', 'ASP', 'ALA', 'TRP', 'LEU', 'SER', 'GLY', 'PRO', 'ASN', 'ALA', 'ASN', 'LYS', 'ILE', 'GLU', 'VAL', 'VAL', 'ILE', 'ALA', 'ASN', 'ASN', 'ASP', 'ALA', 'MET', 'ALA', 'MET', 'GLY', 'ALA', 'VAL', 'GLU', 'ALA', 'LEU', 'LYS', 'ALA', 'HIS', 'ASN', 'LYS', 'SER', 'SER', 'ILE', 'PRO', 'VAL', 'PHE', 'GLY', 'VAL', 'ASP', 'ALA', 'LEU', 'PRO', 'GLU', 'ALA', 'LEU', 'ALA', 'LEU', 'VAL', 'LYS', 'SER', 'GLY', 'ALA', 'LEU', 'ALA', 'GLY', 'THR', 'VAL', 'LEU', 'ASN', 'ASP', 'ALA', 'ASN', 'ASN', 'GLN', 'ALA', 'LYS', 'ALA', 'THR', 'PHE', 'ASP', 'LEU', 'ALA', 'LYS', 'ASN', 'LEU', 'ALA', 'ASP', 'GLY', 'LYS', 'GLY', 'ALA', 'ALA', 'ASP', 'GLY', 'THR', 'ASN', 'TRP', 'LYS', 'ILE', 'ASP', 'ASN', 'LYS', 'VAL', 'VAL', 'ARG', 'VAL', 'PRO', 'TYR', 'VAL', 'GLY', 'VAL', 'ASP', 'LYS', 'ASP', 'ASN', 'LEU', 'ALA', 'GLU', 'PHE']), 1: (1, 'HETATM', 'A', 310, ['BGC']),
            2: (2, 'HETATM', 'A', 311, ['CA']),
            3: (3, 'HETATM', 'A', 312, ['ACT']),
            4: (4, 'HETATM', 'A', 313, ['CO2']),
            5: (5, 'HETATM', 'A', 314, ['GOL'])
        }
    """
    AtomGroup = mda.Universe(pdbfile).select_atoms('all')
    
    GroupList = []
    index = 0
    for residue in AtomGroup.residues:
        if hasattr(residue.atoms, 'record_types'):
            RecordTypes = np.unique(residue.atoms.record_types)
        else:
            RecordTypes = ['ATOM']
        resid = residue.resid
        resname = residue.resname
        ChainId = residue.segid
        ResnameDictionary = pd.read_pickle(__ResidueTypeFile)
        if resname in ResnameDictionary:
            ResType = ResnameDictionary[resname]
        else:
            ResType = 'Other'
        if len(RecordTypes) > 1:
            logging.warning('Found more than one record type in %s-%d.' % (resname, resid))
        if len(GroupList) == 0:
            GroupList = [(index, RecordTypes[0], ChainId, [resid], [resname], ResType)]
        else:
            LastRecordType = GroupList[-1][1]
            LastChainId = GroupList[-1][2]
            LastResid = GroupList[-1][3][-1]
            LastResname = GroupList[-1][4][-1]
            LastResType = GroupList[-1][5]
            if RecordTypes[0] != LastRecordType or ChainId != LastChainId    \
              or (RecordTypes[0] == 'HETATM' and resname != LastResname)    \
              or (resid-LastResid) != 1:# or ResType != LastResType           :
                index += 1
                GroupList.append((index, RecordTypes[0], ChainId, [resid], [resname], ResType))
            else:
                GroupList[-1][3].append(resid)
                GroupList[-1][4].append(resname)
    GroupDict = {}
    for group in GroupList:
        index = group[0]
        RecordType = group[1]
        ChainId = group[2]
        ResidRange = str(group[3][0]) if len(group[3])==1 else "%d:%d"%(group[3][0], group[3][-1])
        resnames = group[4]
        ResType = group[5]
        GroupDict[index] = (index, RecordType, ChainId, ResidRange, resnames, ResType)

    return GroupDict

def print_group(GroupDict:dict) -> list:
        """Print the group information with format.

        Returns:
            list: [List of lines with every line contains a group information]
            example:
               Index   Chain   Record Type     Group Type      Residue Range      Resname
                   0       A          ATOM        Protein              2:306   amino acid
                   1       A        HETATM          Other                310          BGC
                   2       A        HETATM          Other                311           CA
                   3       A        HETATM          Other                312          ACT
                   4       A        HETATM          Other                313          CO2
                   5       A        HETATM          Other                314          GOL
                   6       A        HETATM          Water            501:918          HOH
        """
        title = "Index\tChain\tRecord Type\tGroup Type\tResidue Range\tResname\n"
        LineFormat = "%5d\t%5s\t%11s\t%10s\t%13s\t%10s\n"
        lines = title
        for k,v in GroupDict.items():
            index, RecordType, ChainId, ResidRange, resnames, GroupType = v
            if len(resnames)>2:
                resname = GroupType
            else:
                resname = resnames[0]
            line = LineFormat % (index+1, ChainId, RecordType, GroupType, ResidRange, resname)
            lines += line
        print(lines)

def select_group(pdbfile: str, selectedGroup: list, outfile: str, trajfile=None) -> None:
    """[summary]

    Args:
        pdbfile (str): [description]
        selectedGroup (list): [description]
        outfile (str): [description]
    """

    universe = mda.Universe(pdbfile)
    if trajfile:
        universe.load_new(trajfile)
    groupDict = detect_group(pdbfile)
    groupAtoms = []
    for groupIndex in selectedGroup:
        if groupIndex not in groupDict.keys():
            raise Exception('Not found goup %d in the group list.'%groupIndex, groupDict)
        index, RecordType, ChainId, ResidRange, resnames, ResType = groupDict[groupIndex]
        if RecordType == 'HETATM':
            selection = 'segid %s and resname %s'%(ChainId, resnames[0])
        else:
            if ':' in ResidRange:
                startResid, endResid = ResidRange.split(':')
            else:
                startResid, endResid = ResidRange, ResidRange
            if ChainId.strip():
                selection = 'segid %s and resid %s-%s'%(ChainId, startResid, endResid)
            else:
                selection = 'resid %s-%s'%(startResid, endResid)
        print(selection)
        atoms = universe.select_atoms(selection)
        groupAtoms.append(atoms)
    with mda.Writer(outfile, multiframe=False) as fw:
        for atoms in groupAtoms:
            fw.write(atoms)


def gmx_index_group(pdbfile:str, selectedGroup: list, indexFile='index.ndx') -> None:
    universe = mda.Universe(pdbfile)
    groupDict = detect_group(pdbfile)
    groupAtoms = []
    print(groupDict)
    for groupIndex in selectedGroup:
        if groupIndex not in groupDict.keys():
            raise Exception('Not found goup %d in the group list.'%groupIndex, groupDict)
        index, RecordType, ChainId, ResidRange, resnames, ResType = groupDict[groupIndex]
        if RecordType == 'HETATM':
            selection = 'segid %s and resname %s'%(ChainId, resnames[0])
        else:
            if ':' in ResidRange:
                startResid, endResid = ResidRange.split(':')
            else:
                startResid, endResid = ResidRange, ResidRange
            if ChainId.strip():
                selection = 'segid %s and resid %s-%s'%(ChainId, startResid, endResid)
            else:
                selection = 'resid %s-%s'%(startResid, endResid)
        atoms = universe.select_atoms(selection)
        groupAtoms.append(atoms)

    with mda.selections.gromacs.SelectionWriter(indexFile, mode='w') as ndx:
        receptor, ligand = groupAtoms[0], groupAtoms[1]
        ndx.write(receptor+ligand, name='Complex')
        ndx.write(receptor, name='receptor')
        ndx.write(ligand, name='ligand')


def main():
    parser = argparse.ArgumentParser(description='Detect group information from input file.')
    parser.add_argument('-i', dest='inp', help='Input file with pdb tpr', required=True)

    args = parser.parse_args()

    inp = args.inp
    GroupDict = detect_group(inp)
    print_group(GroupDict)

if __name__ == "__main__":
    main()

