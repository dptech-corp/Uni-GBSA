Uni-GBSA: An Automatic Workflow to Perform MM/GB(PB)SA Calculations for Virtual Screening
==============================================================================

[[ChemRxiv](https://chemrxiv.org/engage/chemrxiv/article-details/63345399f764e656800664e7)]

## Backgroud

Molecular mechanics/Generalized-Born (Poisson–Boltzmann) surface area (MM/GB(PB)SA), which balance accuracy and efficiency, is a good choice for evaluating binding free energy in virtual screening. Uni-GBSA, an automatic workflow to perform MM/GB(PB)SA calculations from force field building, structure optimization to free energy calculation. 

## Install
### Install by conda
To run uni-GBSA, you need install some third-party software( acpype, gmx_MMPBSA, lickit, etc.).
```Bash
conda create -n gbsa -c conda-forge acpype openmpi mpi4py gromacs
conda activate gbsa
pip install unigbsa gmx_MMPBSA>=1.5.6 lickit
```

### Install by dokcer images
You can build a dokcer image by this file or just pull from the docker hub `docker pull dockerymh/unigbsa`
```Plaintext
FROM continuumio/miniconda3

# 1. create a enverioment
SHELL ["/bin/bash", "-c"]
RUN conda create -n gbsa -c conda-forge acpype openmpi mpi4py gromacs \
&&  echo 'conda activate gbsa' >> ~/.bashrc \
&&  rm -rf /opt/conda/pkgs/* 

# 2. install unigbsa
RUN source ~/.bashrc \ 
&&  pip install unigbsa gmx_MMPBSA>=1.5.6 lickit \
&&  rm -rf ~/.cache/*

```

## Usage & Examples

### Usage
```bash
$ unigbsa-pipeline -h
usage: unigbsa-pipeline [-h] -i RECEPTOR [-l LIGAND [LIGAND ...]] [-c CONFIG] [-d LIGDIR] [-f PBSAFILE] [-o OUTFILE] [-nt THREAD]
                        [--decomp] [--verbose] [-v]

GBSA Calculation. Version: 0.0.9_dev

options:
  -h, --help            show this help message and exit
  -i RECEPTOR           Input protein file with pdb format.
  -l LIGAND [LIGAND ...]
                        Ligand files to calculate binding energy.
  -c CONFIG             Configue file, default: /opt/miniconda3/envs/test/lib/python3.10/site-packages/unigbsa-0.0.9.dev0-py3.10.egg/unigbsa/data/default.ini
  -d LIGDIR             Floder contains many ligand files. file format: .mol or .sdf
  -f PBSAFILE           gmx_MMPBSA input file. default=None
  -o OUTFILE            Output file.
  -nt THREAD            Set number of thread to run this program.
  --decomp              Decompose the free energy. default:False
  --verbose             Keep all the files.
  -v, --version         show program's version number and exit
```

### Example
```bash
$ unigbsa-pipeline -i example/1ceb/1ceb_protein.pdb -l example/1ceb/1ceb_ligand.sdf -o BindingEnergy.csv

10/08/2022 13:46:09 PM - INFO - Build protein topology.
10/08/2022 13:46:10 PM - INFO - Build ligand topology: 1ceb_ligand
1 molecule converted
10/08/2022 13:46:13 PM - INFO - Running energy minimization: 1ceb_ligand
10/08/2022 13:46:14 PM - INFO - Run the MMPB(GB)SA.
10/08/2022 13:46:18 PM - INFO - Clean the results.
================================================================================
Results: Energy.csv Dec.csv
Frames    mode    detal_G(kcal/mole)
     1      gb              -20.4421

```

## Other Tools
This packge contains many command lines: `unigbsa-pipeline`, `unigbsa-traj`, `unigbsa-pbc`, `unigbsa-buildtop`, `unigbsa-buildsys`, `unigbsa-md`.

### unigbsa-pipeline
>A very simple pipeline to calculate the PBSA/GBSA value. You just need input a protein file and some ligands files. It will obtain the PBSA/GBSA value for this ligands.


* If you want do minimization or MD simulation for the complex. Just use the ``unigbsa-pipeline``
```Bash
usage: unigbsa-pipeline [-h] -i RECEPTOR [-l LIGAND [LIGAND ...]] [-c CONFIG] [-d LIGDIR] [-f PBSAFILE] [-o OUTFILE] [-nt THREAD] [--decomp] [--verbose]

GBSA Calculation.

optional arguments:
  -h, --help            show this help message and exit
  -i RECEPTOR           Input protein file with pdb format.
  -l LIGAND [LIGAND ...]
                        Ligand files to calculate binding energy.
  -c CONFIG             Configue file, default: default.ini
  -d LIGDIR             Floder contains many ligand files. file format: .mol or .sdf
  -f PBSAFILE           gmx_MMPBSA input file. default=None
  -o OUTFILE            Output file.
  -nt THREAD            Set number of thread to run this program.
  --decomp              Decompose the free energy. default:False
  --verbose             Keep all the files.
```

You can change the parameters for calculations by settig a configue file(`default.ini`). 
```
; parameters for simulation
[simulation]
; input pose process method: 
;   input   -   just use input pose to calculation
;   em      -   run a simple energy minimizaion for the input poses
;   md      -   run a md simulation for the input poses
mode = em

; simulation box type: triclinic, cubic, dodecahedron, octahedron
boxtype = triclinic

; Distance between the solute and the simulation box
boxsize = 0.9

; Specify salt concentration (mol/liter). This will add sufficient ions to reach up to the specified concentration
conc = 0.15

; number of md simulation steps
nsteps = 500000

; number of equilibrium simulation(nvt, npt) steps
eqsteps = 50000

; number of structure to save for the md simulation
nframe = 100

; protein forcefield (gromacs engine)
proteinforcefield = amber03

; ligand forcefield (acpype engine)
ligandforcefield = gaff
; ligand charge method: bcc, gas
ligandCharge = bcc


; parameters for PBSA/GBSA calculation, support all the gmx_MMPBSA parameters
[GBSA]
; calculation name
sys_name = GBSA

; calculation mode, Separated by commas. gb,pb,decomposition
modes = gb

; best parameters for PBSA/GBSA calculations obtained from Wang, Ercheng, et al. Chemical reviews 119.16 (2019): 9478-9508.
igb = 2
indi = 4.0
exdi = 80.0
```

### unigbsa-traj
>Calculate the PBSA/GBSA value for a md trajectory. Most important, you need to prepare a gromacs `index.ndx` file which contains two groups named `RECEPTOR` and `LIGAND`.

````
unigbsa-traj -h
usage: unigbsa-traj [-h] -i INP -p TOP -ndx NDX [-m {gb,pb,pb+gb,gb+pb}] [-t TRAJ] [-indi INDI] [-dec] [-D]

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

### unigbsa-buildtop
>Build topology for protein or ligand by gromacs.
```Bash
unigbsa-buildtop -h
usage: unigbsa-buildtop [-h] [-p PROTEIN] [-l LIGAND] [-pf PROTFORCE] [-lf {gaff,gaff2}] [-o OUTDIR] [-c] [-verbose]

Build topology file for input file.

optional arguments:
  -h, --help        show this help message and exit
  -p PROTEIN        Protein file or directory to build topology.
  -l LIGAND         Ligand file or directory to build topology.
  -pf PROTFORCE     Protein forcefield.
  -lf {gaff,gaff2}  Ligand forcefiled: gaff or gaff2.
  -o OUTDIR         A output directory.
  -c                Combine the protein and ligand topology. Suppport for one protein and more ligands. default:True
  -verbose          Keep the directory or not.
```
### unigbsa-buildsys
>Build simulation system for protein or ligand.
```bash
unigbsa-buildsys -h
usage: unigbsa-buildsys [-h] -p PROTEIN [-l LIGAND] [-pf PROTFORCE] [-lf {gaff,gaff2}] [-bt BOXTYPE] [-box BOX BOX BOX] [-d D] [-conc CONC] [-o OUTDIR]

Build MD simulation for input file.

optional arguments:
  -h, --help        show this help message and exit
  -p PROTEIN        Protein file for the simulation.
  -l LIGAND         Ligand file or directory for the simulation.
  -pf PROTFORCE     Protein forcefield.
  -lf {gaff,gaff2}  Ligand forcefiled: gaff or gaff2.
  -bt BOXTYPE       Simulation box type, default: triclinic
  -box BOX BOX BOX  Simulation box size.
  -d D              Distance between the solute and the box.
  -conc CONC        Specify salt concentration (mol/liter). default=0.15
  -o OUTDIR         A output directory.
```

### unigbsa-md
>Run a MD simulation for input protein or ligand.
```Bash
unigbsa-md -h
usage: unigbsa-md [-h] -p PROTEIN [-l LIGAND] [-pf PROTFORCE] [-lf {gaff,gaff2}] [-bt BOXTYPE] [-box BOX BOX BOX] [-d D] [-conc CONC] [-o OUTDIR] [-nstep NSTEP] [-nframe NFRAME] [-verbose]

Run MD simulation for input file.

optional arguments:
  -h, --help        show this help message and exit
  -p PROTEIN        Protein file for the simulation.
  -l LIGAND         Ligand file or directory for the simulation.
  -pf PROTFORCE     Protein forcefield.
  -lf {gaff,gaff2}  Ligand forcefiled: gaff or gaff2.
  -bt BOXTYPE       Simulation box type, default: triclinic
  -box BOX BOX BOX  Simulation box size.
  -d D              Distance between the solute and the box.
  -conc CONC        Specify salt concentration (mol/liter). default=0.15
  -o OUTDIR         A output directory.
  -nstep NSTEP      Simulation steps. default:2500
  -nframe NFRAME    Number of frame to save for the xtc file. default:100
  -verbose          Keep all the files in the simulation.
```

### unigbsa-pbc
>Process PBC condition for the gromacs trajectory.
```Bash
unigbsa-pbc -h
usage: unigbsa-pbc [-h] -s TPR -f XTC [-o OUT] [-n NDX]

Auto process PBC for gromacs MD trajector.

optional arguments:
  -h, --help  show this help message and exit
  -s TPR      TPR file generated from gromacs or coordinate file.
  -f XTC      Trajector file to process PBC.
  -o OUT      Result file after processed PBC.
  -n NDX      Index file contains the center and output group.
```


### More Examples

* Give a protein and some ligand files. Obtain the PBSA with ``unigbsa-pipeline``
````Bash
unigbsa-pipeline -i ./example/2fvy/protein.pdb -l ./example/2fvy/BGC.mol2
````

* Calculate PBSA value with ``unigbsa-traj``
```Bash
unigbsa-traj -i example/3f/complex.pdb -p example/3f/complex.top -ndx example/3f/index.ndx -m pb gb -t example/3f/complex.pdb
```

* Build topology for protein or ligand by gromacs. ``unigbsa-buildtop``
```bash
unigbsa-buildtop -p example/2fvy/protein.pdb -pf amber99sb -o topol  # build gromacs topology for a single protein
unigbsa-buildtop -p example/2fvy/protein.pdb -pf amber99sb -l example/2fvy/BGC.mol2 -lf gaff -o 2fvy_topol -c # build gromacs topology for protein and ligand complex
```

* Build simulation system with ``unigbsa-buildsys``

* Run MD simulation with ``unigbsa-md``

* Process PBC condition with ``unigbsa-pbc``


## Citation
```plaintext
Yang M, Wang D, Zheng H. Uni-GBSA: An Automatic Workflow to Perform MM/GB(PB)SA Calculations for Virtual Screening. ChemRxiv. Cambridge: Cambridge Open Engage; 2022; This content is a preprint and has not been peer-reviewed.
```