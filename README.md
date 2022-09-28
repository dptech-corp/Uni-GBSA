# Uni-GBSA
## Backgroud

Molecular mechanics/Generalized-Born (Poisson–Boltzmann) surface area (MM/GB(PB)SA), which balance accuracy and efficiency, is a good choice for evaluating binding free energy in virtual screening. Uni-GBSA, an automatic workflow to perform MM/GB(PB)SA calculations from force field building, structure optimization to free energy calculation. 

## Install
### Install by conda
```Bash
conda create -n amber21 -c conda-forge ambertools=21 acpype=2021.02 openmpi mpi4py gromacs
conda activate amber21
pip install gmx_MMPBSA
git clone https://github.com/dptech-corp/Uni-GBSA.git
cd Uni-GBSA
python setup.py install
```

### Install by dokcer images
You can build a dokcer image by this file or just pull from the docker hub `docker pull dockerymh/hmtpbsa:gmx_openmpi`
```Plaintext
FROM continuumio/miniconda3

# 1. create a enverioment
SHELL ["/bin/bash", "-c"]
RUN conda create -n amber21 -c conda-forge ambertools=21 acpype=2021.02 gromacs openmpi mpi4py \
&&  echo 'conda activate amber21' >> ~/.bashrc \
&&  rm -rf /opt/conda/pkgs/* 


# 2. install requirments
RUN source ~/.bashrc \ 
&&  pip install gmx_MMPBSA \
&&  rm -rf ~/.cache/*

# 3. install hmtpbsa
RUN git clone https://github.com/dptech-corp/Hermite-MMPBSA.git \
&&  cd Hermite-MMPBSA \
&&  python setup.py install \
&&  cd ../ \
&&  rm -rf Hermite-MMPBSA
```

## Usage & Examples
This packge contains many command lines: `hmtpbsa-pipeline`, `hmtpbsa-traj`, `hmtpbsa-pbc`, `hmtpbsa-buildtop`, `hmtpbsa-buildsys`, `hmtpbsa-md`.

### hmtpbsa-pipeline
>A very simple pipeline to calculate the PBSA/GBSA value. You just need input a protein file and some ligands files. It will obtain the PBSA/GBSA value for this ligands.


* If you want do minimization or MD simulation for the complex. Just use the ``hmtpbsa-pipeline``
```Bash
usage: hmtpbsa-pipeline [-h] -i RECEPTOR [-l LIGAND [LIGAND ...]] [-c CONFIG] [-d LIGDIR] [-f PBSAFILE] [-o OUTFILE] [-nt THREAD] [--decomp] [--verbose]

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
[PBSA]
; calculation name
sys_name = GBSA

; calculation mode, Separated by commas. gb,pb,decomposition
modes = gb

; best parameters for PBSA/GBSA calculations obtained from Wang, Ercheng, et al. Chemical reviews 119.16 (2019): 9478-9508.
igb = 2
indi = 4.0
exdi = 80.0
```

### hmtpbsa-traj
>Calculate the PBSA/GBSA value for a md trajectory. Most important, you need to prepare a gromacs `index.ndx` file which contains two groups named `RECEPTOR` and `LIGAND`.

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

### hmtpbsa-buildtop
>Build topology for protein or ligand by gromacs.
```Bash
hmtpbsa-buildtop -h
usage: hmtpbsa-buildtop [-h] [-p PROTEIN] [-l LIGAND] [-pf PROTFORCE] [-lf {gaff,gaff2}] [-o OUTDIR] [-c] [-verbose]

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
### hmtpbsa-buildsys
>Build simulation system for protein or ligand.
```bash
hmtpbsa-buildsys -h
usage: hmtpbsa-buildsys [-h] -p PROTEIN [-l LIGAND] [-pf PROTFORCE] [-lf {gaff,gaff2}] [-bt BOXTYPE] [-box BOX BOX BOX] [-d D] [-conc CONC] [-o OUTDIR]

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

### hmtpbsa-md
>Run a MD simulation for input protein or ligand.
```Bash
hmtpbsa-md -h
usage: hmtpbsa-md [-h] -p PROTEIN [-l LIGAND] [-pf PROTFORCE] [-lf {gaff,gaff2}] [-bt BOXTYPE] [-box BOX BOX BOX] [-d D] [-conc CONC] [-o OUTDIR] [-nstep NSTEP] [-nframe NFRAME] [-verbose]

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

### hmtpbsa-pbc
>Process PBC condition for the gromacs trajectory.
```Bash
hmtpbsa-pbc -h
usage: hmtpbsa-pbc [-h] -s TPR -f XTC [-o OUT] [-n NDX]

Auto process PBC for gromacs MD trajector.

optional arguments:
  -h, --help  show this help message and exit
  -s TPR      TPR file generated from gromacs or coordinate file.
  -f XTC      Trajector file to process PBC.
  -o OUT      Result file after processed PBC.
  -n NDX      Index file contains the center and output group.
```


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

* Build simulation system with ``hmtpbsa-buildsys``

* Run MD simulation with ``hmtpbsa-md``

* Process PBC condition with ``hmtpbsa-pbc``
