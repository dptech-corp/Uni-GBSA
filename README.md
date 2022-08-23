# Hermite-MMPBSA
## Backgroud

see the [MM/GB(PB)SA介绍](!https://dptechnology.feishu.cn/wiki/wikcnfUDQ1sL2oXAl5GVDQhzzSb)


## Install
### step 1: install gromacs
```Bash
apt-get update
apt install gromacs
```

### step 2: install ambertools
```Bash
conda create -n amber21 -c conda-forge ambertools=21 acpype=2021.02
conda activate amber21
pip install gmx_MMPBSA==1.5.2

```

### step 3: install hermite-mmpbsa
```Bash
python setup.py install
```

## Usage

The `OMP_NUM_THREADS` environment variable specifies the number of threads to use for parallel regions. If you do not set `OMP_NUM_THREADS` , the number of processors available is the default value `4`


* If you want do minimization or MD simulation for the complex. Just use the ``hmtpbsa-pipeline``
```Bash
hmtpbsa-pipeline -h
usage: hmtpbsa-pipeline [-h] -i RECEPTOR [-l LIGAND [LIGAND ...]] [-c CONFIG] [-d LIGDIR] [-f PBSAFILE] [-o OUTFILE]

GBSA Calculation.

optional arguments:
  -h, --help            show this help message and exit
  -i RECEPTOR           Input protein file with pdb format.
  -l LIGAND [LIGAND ...]
                        Ligand files to calculate binding energy.
  -c CONFIG             Configue file, default: /opt/anaconda3/envs/amber/lib/python3.8/site-packages/hmtpbsa-0.0.2-py3.8.egg/hmtpbsa/data/default.ini
  -d LIGDIR             Floder contains many ligand files. file format: .mol or .sdf
  -f PBSAFILE           gmx_MMPBSA input file. default=None
  -o OUTFILE            Output file.
```

* If you have the gromacs topology and index files. Just use the ``hmtpbsa-traj``
````
hmtpbsa-traj -h
usage: hmtpbsa-traj [-h] -i INP -p TOP -ndx NDX [-m {gb,pb,pb+gb,gb+pb}] [-t TRAJ] [-indi INDI] [-dec] [-D]

Free energy calcaulated by PBSA method.

optional arguments:
  -h, --help            show this help message and exit
  -i INP                A pdb file or a tpr file for the trajectory.
  -p TOP                Gromacs topol file for the system.
  -ndx NDX              Gromacs index file, must contain recepror and ligand group.
  -m {gb,pb,pb+gb,gb+pb}
                        Method to calculate: gb, pb, pb+gb. default:gb
  -t TRAJ               A trajectory file contains many structure frames. File format: xtc, pdb, gro...
  -indi INDI            External dielectric constant. detault: 1.0
  -dec                  Decompose the energy. default:false
  -D                    DEBUG model, keep all the files.
````


## Example

* Give a protein and some ligand files. Obtain the PBSA with ``hmtpbsa-pipeline``
````Bash
hmtpbsa-pipeline -i ./example/2fvy/protein.pdb -l ./example/2fvy/BGC.mol2
````

* Calculate PBSA value with ``hmtpbsa-traj``
```Bash
hmtpbsa-traj -i example/3f/complex.pdb -p example/3f/complex.top -ndx example/3f/index.ndx -m pb gb -t example/3f/complex.pdb
```

* Build topology for protein or ligand by gromacs. ``hmtpbsa-buildtop``
```bash
hmtpbsa-buildtop -p example/2fvy/protein.pdb -pf amber99sb -o topol  # build gromacs topology for a single protein
hmtpbsa-buildtop -p example/2fvy/protein.pdb -pf amber99sb -l example/2fvy/BGC.mol2 -lf gaff -o 2fvy_topol -c # build gromacs topology for protein and ligand complex
```

* Run MD simulation with ``hmtpbsa-md``

* Process PBC condition with ``hmtpbsa-pbc``

* Build simulation system with ``hmtpbsa-buildsys``


* Build topology for protein or ligand by gromacs. ``hmtpbsa-buildtop``
```bash
hmtpbsa-buildtop -p example/2fvy/protein.pdb -pf amber99sb -o topol  # build gromacs topology for a single protein
hmtpbsa-buildtop -p example/2fvy/protein.pdb -pf amber99sb -l example/2fvy/BGC.mol2 -lf gaff -o 2fvy_topol -c # build gromacs topology for protein and ligand complex
```

* Run MD simulation with ``hmtpbsa-md``

* Process PBC condition with ``hmtpbsa-pbc``

* Build simulation system with ``hmtpbsa-buildsys``
