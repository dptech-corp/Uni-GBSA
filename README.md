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
```Bash
hmtmmpbsa -h
Usage: hmtmmpbsa [-h] -i INP [-t TRAJ] [-o OUTP] [-ndx NDX] [-D]

Free energy calcaulated by MMPBSA method.

optional arguments:
  -h, --help  show this help message and exit
  -i INP      A pdb file or a tpr file to calculate the free energy.
  -t TRAJ     A trajectory file contains many structure frames. File format: xtc, pdb, gro...
  -o OUTP     Output floder to save results.
  -ndx NDX    Index file.
  -D          DEBUG model, keep all the files.
```

### Single PDB file
```Bash
hmtmmpbsa -i example/2fvy.pdb -o single-pdb << EOF
1,2
EOF
```

### Trajectory
```Bash
hmtmmpbsa -i example/ST/com.tpr -t example/ST/com_traj.xtc -o trajs-test << EOF
1,2
EOF

```