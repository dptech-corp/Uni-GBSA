from pdb4amber.pdb4amber import run as amber_prep
from pdb2pqr.main import main_driver
import argparse
def parse_args(**kwargs):
    return kwargs

def preparation(pdbfile: str, outfile: str) -> None:
    ArgsDict = {
        "arg_pdbout":outfile,
        "arg_pdbin":pdbfile,
        "arg_nohyd":True,
        "arg_dry":False,
        "arg_prot":False,
        "arg_amber_compatible_residues":False,
        "arg_strip_atom_mask":None,
        "arg_mutate_string":None,
        "arg_constph":False,
        "arg_mostpop":False,
        "arg_reduce":False,
        "arg_no_reduce_db":False,
        "arg_model":0,
        "arg_add_missing_atoms":True, # A bug when set true
        "arg_elbow":False,
        "arg_logfile":'pdb4amber.log',
        "arg_keep_altlocs":False,
        "arg_leap_template":False,
        "arg_conect":True,
        "arg_noter":False,
    }
    amber_prep(**ArgsDict)

def clean_header(protfile:str):
    '''
    clean the header lines from the pdb file.
    :param protfile: protein file
    :return:
    '''
    headers = []; HF = True
    with open(protfile) as fr:
        lines = fr.readlines()
    with open(protfile, 'w') as fw:
            for line in lines:
                if HF and not (line.startswith('ATOM') or line.startswith('HETATM')):
                    headers.append(line)
                    continue
                else:
                    HF = False
                fw.write(line)
    return headers

def add_header(protfile:str, headers:list):
    '''
    Add header lines for pdb file.
    :param protfile: str, input protein file
    :param headers: list, list line of headers
    :return:
    '''
    if len(headers)==0:
        return 0
    with open(protfile) as fr:
        lines = fr.readlines()
    with open(protfile, 'w') as fw:
        fw.writelines(headers)
        fw.writelines(lines)

def pdb2pqr(pdbfile:str, outfile:str, ffout='AMBER'):
    outpqr = outfile[:-4] + '.pqr'
    argsdict = parse_args(
        alignment=None,
        apbs_input=None, 
        assign_only=False, 
        chains=None, 
        clean=False, 
        debump=True, 
        display_coupled_residues=False, 
        drop_water=False, 
        ff='PARSE', 
        ffout=ffout, 
        filenames=[],
        grid=(0.0, 14.0, 0.1),
        include_header=False, 
        input_path=pdbfile, 
        keep_chain=False, 
        keep_protons=False, 
        ligand=None, 
        log_level='INFO', 
        mutations=None, 
        mutator=None, 
        mutator_options=None, 
        neutralc=False, 
        neutraln=False, 
        opt=True, 
        output_pqr=outpqr, 
        pH=7.0, 
        parameters='/opt/anaconda3/envs/amber/lib/python3.8/site-packages/propka/propka.cfg', 
        pdb_output=outfile, 
        ph=7.0, 
        pka_method=None, 
        protonate_all=False, 
        reference='neutral', 
        reuse_ligand_mol2_file=False, 
        thermophiles=None, 
        titrate_only=None, 
        userff=None, 
        usernames=None, 
        whitespace=False, 
        window=(0.0, 14.0, 1.0)
    )
    headers = clean_header(pdbfile)
    argsdict = argparse.Namespace(**argsdict)
    _ = main_driver(argsdict)
    add_header(pdbfile, headers)

def main():
    parser = argparse.ArgumentParser(description='pdb2pqr30')
    parser.add_argument('-i', dest='inp', help='Input pdb file.', required=True)
    parser.add_argument('-o', dest='oup', help='Output pdb file.', required=True)

    args = parser.parse_args()
    pdbfile, outfile = args.inp, args.oup
    pdb2pqr(pdbfile, outfile)


if __name__ == "__main__":
    main()
