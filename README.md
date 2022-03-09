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
conda create -n amber21 -c conda-forge ambertools=21

```

### step 3: install requirment
```Bash
conda activate amber21
conda install mdanalysis==2.0.0 gmx_MMPBSA pdb2pqr -c conda-forge
```

### step 4: installl hermite-mmpbsa
```Bash
python setup.py install
```

## Usage
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


* If you want do minimization or MD simulation for the complex. Just use the ``hmtpbsa-pipeline``
```Bash
hmtpbsa-pipeline -h
usage: hmtpbsa-pipeline [-h] -i RECEPTOR -l LIGAND [LIGAND ...] [-c CONFIG]

GBSA Calculation.

optional arguments:
  -h, --help            show this help message and exit
  -i RECEPTOR           Input protein file with pdb format.
  -l LIGAND [LIGAND ...]
                        Ligand files to calculate binding energy.
  -c CONFIG             Configue file, default: /opt/anaconda3/envs/amber/lib/python3.8/site-
                        packages/hmtpbsa-0.0.2-py3.8.egg/hmtpbsa/data/detault.ini
```

## Example
* Calculate PBSA value with ``hmtpbsa-traj``
```Bash
hmtpbsa-traj -i example/3f/complex.pdb -p example/3f/complex.top -ndx example/3f/complex.ndx -m pb+gb -t example/3f/complex.pdb
```

* Give a protein and some ligand files. Obtain the binding energy with ``hmtpbsa-pipeline``
````Bash
hmtpbsa-pipeline -i example/md/protein.pdb -l example/md/3f.mol
````
