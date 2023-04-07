from enum import Enum

import sh
import os
import json
from pathlib import Path
from dp.launching.report import Report, AutoReportElement, ReportSection
from dp.launching.cli import SubParser, run_sp_and_exit, default_minimal_exception_handler
from dp.launching.typing.basic import BaseModel, Field, Float, Set, List
from dp.launching.typing.io import InputFilePath, OutputDirectory


modedic = {
    'GB-HCT': ('GB', 1, 'gb-1'),
    'GB-OBC1': ('GB', 2, 'gb-2'),
    'GB-OBC2': ('GB', 5, 'gb-5'),
    'GB-Neck': ('GB', 7, 'gb-7'),
    'GB-Neck2': ('GB', 8, 'gb-8'),
    'PB-1': ('PB', 1, 'pb-1'),
    'PB-2': ('PB', 2, 'pb-2'),
}


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


class SCANUploadFilesModel(BaseModel):
    experiment_info: InputFilePath = \
        Field(ftypes=['csv'], description='Input experiment file.')
    input_protein: InputFilePath = \
        Field(ftypes=['pdb'],
              description="Input protein file with pdb format.")
    input_ligands: List[InputFilePath] = \
        Field(ftypes=['sdf'],
              description='Ligand files to perform the calculation.')


class SCANArgs(BaseModel):
    modes: Set[AlgorithmMode] = Field(default=AlgorithmMode.gb2, render_type="radio",
                                      description="MM PB/GBSA calculation mode.")
    indi: str = Field(default="1.0, 2.0, 2.5, 3.0, 3.5, 4.0",
                      description='Internal dielectric constant')
    endi: Float = Field(default=80.0,
                        description='External dielectric constant')
    protein_forcefield: Set[AlgorithmProteinForcefield] =\
        Field(default=AlgorithmProteinForcefield.amber03,
              description='Protein forcefiled.')
    ligand_forcefield: Set[AlgorithmLigandForcefield] =\
        Field(default=AlgorithmLigandForcefield.gaff,
              description='Ligand forcefield.')
    ligand_charge: Set[AlgorithmLigandCharge] =\
        Field(default=AlgorithmLigandCharge.bcc,
              description='Ligand charge method.')


class SCANDefaultModel(SCANArgs, BaseModel):
    pass


class SCANParamsModel(BaseModel):
    pass


class Output(BaseModel):
    output_dir: OutputDirectory


class SCANModel(SCANUploadFilesModel,
                SCANDefaultModel,
                SCANParamsModel,
                Output,
                BaseModel):
    pass


class GBSAUploadFilesModel(BaseModel):
    input_protein: InputFilePath = \
        Field(ftypes=['pdb'],
              description="Input protein file with pdb format.")
    input_ligands: List[InputFilePath] = \
        Field(ftypes=['sdf'],
              description='Ligand files to perform the calculation.')


class AlgorithmOptions(BaseModel):
    protein_forcefield = Field(default=AlgorithmProteinForcefield.amber03,
                               description='Protein forcefiled.')
    ligand_forcefield = Field(default=AlgorithmLigandForcefield.gaff,
                              description='Ligand forcefield.')
    ligand_charge = Field(default=AlgorithmLigandCharge.bcc,
                          description='Ligand charge method.')
    mode = Field(default=AlgorithmMode.gb2,
                 description="MM PB/GBSA calculation mode.")
    indi: Float = Field(default=4.0,
                        description='Internal dielectric constant')
    endi: Float = Field(default=80.0,
                        description='External dielectric constant')


class GBSAModel(
    Output,
    GBSAUploadFilesModel,
    AlgorithmOptions,
    BaseModel
  ):
    ...


def gbsa_runner(opts: GBSAModel) -> int:
    status = 0
    try:
        from pprint import pprint
        print('Opts:')
        pprint(opts.dict())
        modes = {
            'GB-HCT': ('gb', 1),
            'GB-OBC1': ('gb', 2),
            'GB-OBC2': ('gb', 5),
            'GB-Neck': ('gb', 7),
            'GB-Neck2': ('gb', 8),
            'PB-1': ('pb', 1),
            'PB-2': ('pb', 2),
        }
        default = {
            'simulation': {
                'mode': 'em',
                'boxtype': 'triclinic',
                'boxsize': 0.9,
                'conc': 0.15,
                'nsteps': 500000,
                'nframe': 100,
                'eqsteps': 50000,
                'proteinforcefield': opts.protein_forcefield,
                'ligandforcefield': opts.ligand_forcefield,
                'maxsol': 0,
                'ligandCharge': opts.ligand_charge
                },
            'GBSA': {
                'sys_name': 'GBSA',
                'modes': modes[opts.mode][0],
                'igb': modes[opts.mode][1],
                'indi': opts.indi,
                'exdi': opts.endi
                }
            }
        input_protein = opts.input_protein.get_path()
        input_dir = os.path.dirname(input_protein)
        output_dir = str(opts.output_dir)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        csvoutfile = os.path.join(output_dir, 'GBSA-output.csv')
        configfile = Path(os.path.join(input_dir, 'config.json'))
        with open(configfile,  'w') as fw:
            json.dump(default, fw)

        cmd = sh.Command('unigbsa-pipeline')
        args = [
            '-i', opts.input_protein.get_path(),
            '-l', " ".join(opts.input_ligands),
            '-c', configfile,
            '-o', csvoutfile
        ]

        for line in cmd(*args, _iter=True):
            print(line)
        with open(csvoutfile) as fr:
            lines = fr.readlines()
        with open(csvoutfile, 'w') as fw:
            for line in lines:
                llist = line.strip().split(',')
                newlist = [llist[0]]
                for li in llist[2:-1]:
                    newlist.append(str(round(float(li), 4)))
                newline = ','.join(newlist) + '\n'
                fw.write(newline)
        table_section = ReportSection(
                title='GBSA output',
                elements=[
                    AutoReportElement(title="",
                                      path='GBSA-output.csv',
                                      description='')
                ]
            )
        report = Report(title='Results', sections=[table_section],
                        description='')
        report.save(output_dir)
    except Exception as e:
        print(e)
        status = -1
    print(f"job status: {status}")
    return status


def scan_runner(opts: GBSAModel) -> int:
    status = 0
    try:
        from pprint import pprint
        print('Opts:')
        pprint(opts.dict())
        default = {
            "simulation": {
                "mode": "em",
                "boxtype": "triclinic",
                "boxsize": 0.9,
                "conc": 0.15,
                "nsteps": 500000,
                "nframe": 100,
                "eqsteps": 50000,
                "proteinforcefield": [pf.value for pf in opts.protein_forcefield],
                "ligandforcefield": [lf.value for lf in opts.ligand_forcefield],
                "maxsol": 0,
                "ligandCharge": [lc.value for lc in opts.ligand_charge]
            },
            "GBSA": {
                "sys_name": "GBSA",
                "modes": [modedic[k][2] for k in opts.modes],
                "indi": [float(indi) for indi in opts.indi.split(',')],
                "exdi": opts.endi,
                "nonpolarsurfConst": "0.0",
                "surften": "0.0072"
            }
        }

        input_protein = opts.input_protein.get_path()
        input_dir = os.path.dirname(input_protein)
        configfile = Path(os.path.join(input_dir, 'config.json'))
        with open(configfile,  'w') as fw:
            json.dump(default, fw, indent='  ')
        cmd = sh.Command('unigbsa-scan')
        args = [
            '-i', input_protein,
            '-l', opts.input_ligands,
            '-e', opts.experiment_info.get_path(),
            '-c', configfile,
            '-o', opts.output_dir,
        ]
        for line in cmd(*args, _iter=True):
            print(line)
        csvfile = os.path.join(opts.output_dir, 'paras_performance.csv')
        with open(csvfile) as fr:
            lines = fr.readlines()
        with open(csvfile, 'w') as fw:
            for i, line in enumerate(lines):
                llist = line.strip().split(',')
                if i == 1:
                    best_parafile = llist[-1]
                else:
                    newlist = [llist[0]]
                    for li in llist[1:-1]:
                        newlist.append(str(round(float(li), 3)))
                newline = ','.join(newlist) + '\n'
                fw.write(newline)
        with open(best_parafile) as fr:
            js = json.load(fr)
        outfile = os.path.join(opts.output_dir, 'best_paras.json')
        outjs = {
            'simulation': js['simulation'],
            'GBSA': js['GBSA'],
            'results': js['results']
            }
        with open(outfile, 'w') as fw:
            json.dump(outjs, fw)
        table_section = ReportSection(
                title='Scan output',
                elements=[
                    AutoReportElement(title="",
                                      path='paras_performance.csv',
                                      description='')
                ]
            )
        json_section = ReportSection(
            title='Best parameters',
            elements=[AutoReportElement(title='',
                                        path='best_paras.json',
                                        description='')]
        )
        report = Report(title='Results', sections=[table_section, json_section],
                        description='')
        output_dir = str(opts.output_dir)
        report.save(output_dir)
    except Exception as e:
        print(e)
        status = -1
    print(f"job status: {status}")
    return status


def to_parser():
    from collections import OrderedDict
    return OrderedDict({
        "gbsa": SubParser(GBSAModel, gbsa_runner, "Run GBSA Model"),
        "scan": SubParser(SCANModel, scan_runner, "Run SCAN Model"),
    })


if __name__ == '__main__':
    print("current file __name__ == __main__")

    # to_parser()(sys.argv[1:])
    run_sp_and_exit(to_parser(), description="Uni-GBSA", version="0.1.3",
                    exception_handler=default_minimal_exception_handler)
