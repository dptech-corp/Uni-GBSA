import sh
import json
import uuid
from pathlib import Path
from dp.launching.typing import BaseModel, Field

from dp.launching.typing import InputFilePath, OutputDirectory
from dp.launching.typing import List, Enum, Float, Set

from dp.launching.cli import to_runner, default_minimal_exception_handler


class IOOptions(BaseModel):
    input_protein: InputFilePath = Field(..., ftypes=['pdb'],
                                         description="Input protein file with pdb format.")
    input_ligands: List[InputFilePath] = Field(..., ftypes=['sdf'],
                                               description='Ligand files to perform the calculation.')
    input_exp: InputFilePath = Field(..., ftypes=['csv'],
                                     description='Experiment Info.')
    output_dir: OutputDirectory = Field(default="./output")


class AlgorithmProteinForcefield(str, Enum):
    """
    Protein forcefiled enumerate
    """

    amber03 = 'amber03'
    amber99sb = "amber99sb"
    amber99sb_ildn = 'amber99sb-ildn'


class AlgorithmLigandForcefield(str, Enum):
    """
    Ligand forcefield enumerate
    """
    gaff = "gaff"
    gaff2 = "gaff2"


class AlgorithmLigandCharge(str, Enum):
    """
    Ligand charge method
    """
    bcc = "bcc"
    gaseiger = "gasteiger"
    mmff94 = "mmff94"


class AlgorithmMode(str, Enum):
    """
    MM PB/GBSA mode
    """
    gb1 = 'GB-HCT'
    gb2 = 'GB-OBC1'
    gb5 = 'GB-OBC2'
    gb7 = 'GB-Neck'
    gb8 = 'GB-Neck2'
    pb1 = 'PB-1'
    pb2 = 'PB-2'


class AlgorithmOptions(BaseModel):
    protein_forcefield: Set(AlgorithmLigandForcefield) =\
        Field(default=AlgorithmProteinForcefield.amber03,
              description='Protein forcefiled.')
    ligand_forcefield: Set(AlgorithmLigandForcefield) =\
        Field(default=AlgorithmLigandForcefield.gaff,
              description='Ligand forcefield.')
    ligand_charge: Set(AlgorithmLigandCharge) =\
        Field(default=AlgorithmLigandCharge.bcc,
              description='Ligand charge method.')
    mode: Set(AlgorithmMode) = \
        Field(default=AlgorithmMode.gb2,
              description="MM PB/GBSA calculation mode.")
    indi: Float = Field(default=4.0,
                        description='Internal dielectric constant')
    endi: Float = Field(default=80.0,
                        description='External dielectric constant')


class GlobalOptions(IOOptions, AlgorithmOptions, BaseModel):
    ...


def pipeline_runner(opts: GlobalOptions) -> int:
    status = 0
    from pprint import pprint
    print('Opts:')
    pprint(opts.dict())
    modes = {
        'GB-HCT': ('GB', 1),
        'GB-OBC1': ('GB', 2),
        'GB-OBC2': ('GB', 5),
        'GB-Neck': ('GB', 7),
        'GB-Neck2': ('GB', 8),
        'PB-1': ('PB', 1),
        'PB-2': ('PB', 2),
    }
    default = {
        "simulation": {
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
        },
        "GBSA": {
            "sys_name": "GBSA",
            "modes": ["gb-1", "gb-2", "gb-5", "pb-1", "pb-2"],
            "indi": [2.0, 4.0, 3.0],
            "exdi": "80.0",
            "nonpolarsurfConst": "0.0",
            "surften": "0.0072"
        }
    }
    configfile = Path('/tmp/' + str(uuid.uuid4()) + '.json')
    with open(configfile,  'w') as fw:
        json.dump(default, fw)
    cmd = sh.Command('/opt/conda/envs/amber21/bin/unigbsa-scan')
    args = [
        '-i', opts.input_protein.get_path(),
        '-l', " ".join(opts.input_ligands),
        '-e', opts.input_exp.get_path(),
        '-c', configfile,
        '-o', opts.output_dir,
    ]
    for line in cmd(*args, _iter=True):
        print(line)
    return status


__all__ = ['GlobalOptions']


def to_parser():
    return to_runner(
        GlobalOptions,
        pipeline_runner,
        version='0.1.0',
        exception_handler=default_minimal_exception_handler,
    )


if __name__ == '__main__':
    import sys

    to_parser()(sys.argv[1:])
