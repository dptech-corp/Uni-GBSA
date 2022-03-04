
import os
import sys
import shlex
import logging
import argparse
import numpy as np
from os.path import split

import MDAnalysis as mda

try:
    from .lib.prepare import pdb2pqr
    from .lib.parameters import generate_input_file
    from .lib.groups import detect_group, gmx_index_group, print_group, select_group
except:
    from lib.prepare import pdb2pqr
    from lib.parameters import generate_input_file
    from lib.groups import detect_group, gmx_index_group, print_group, select_group

try:
    from GMXMMPBSA.exceptions import InputError, CommandlineError, GMXMMPBSA_WARNING
    from GMXMMPBSA.infofile import InfoFile
    from GMXMMPBSA import main as pbsaMain
    from GMXMMPBSA.commandlineparser import anaparser, testparser
except ImportError:
    import os
    amberhome = os.getenv('AMBERHOME') or '$AMBERHOME'
    raise ImportError('Could not import Amber Python modules. Please make sure '
                      'you have sourced %s/amber.sh (if you are using sh/ksh/'
                      'bash/zsh) or %s/amber.csh (if you are using csh/tcsh)' %
                      (amberhome, amberhome))

def finalize(app):
    """ We are done. Finish up timers and print out timing info """
    app.timer.done()
    if not app.master:
        app.MPI.Finalize()
        return 0
    app.stdout.write('\nTiming:\n')
    app.timer.print_('setup_gmx', app.stdout)
    app.timer.print_('setup', app.stdout)
    if not app.FILES.rewrite_output:
        app.timer.print_('cpptraj', app.stdout)
        if app.INPUT['alarun']:
            app.timer.print_('muttraj', app.stdout)
        app.timer.print_('calc', app.stdout)
        app.stdout.write('\n')
        if app.INPUT['gbrun']:
            app.timer.print_('gb', app.stdout)
        if app.INPUT['pbrun']:
            app.timer.print_('pb', app.stdout)
        if app.INPUT['nmoderun']:
            app.timer.print_('nmode', app.stdout)
        if app.INPUT['qh_entropy']:
            app.timer.print_('qh', app.stdout)
        app.stdout.write('\n')
    app.timer.print_('output', app.stdout)
    app.timer.print_('global', app.stdout)
    app.remove(app.INPUT['keep_files'])
    logging.info('\n\nThank you for using gmx_MMPBSA. Please cite us if you publish this work with this '
                 'reference:\n    Mario S. Valdés Tresanco, Mario E. Valdes-Tresanco, Pedro A. Valiente, & '
                 'Ernesto Moreno Frías \n    gmx_MMPBSA (Version v1.4.2). '
                 'Zenodo. http://doi.org/10.5281/zenodo.4569307'
                 '\n\nAlso consider citing MMPBSA.py\n    Miller III, B. R., McGee Jr., '
                 'T. D., Swails, J. M. Homeyer, N. Gohlke, H. and Roitberg, A. E.\n    J. Chem. Theory Comput., '
                 '2012, 8 (9) pp 3314-3321\n')
    app.MPI.Finalize()
    end = 0
    if app.FILES.gui:
        import subprocess
        logging.info('Opening gmx_MMPBSA_ana to analyze results...\n')
        g = subprocess.Popen(['gmx_MMPBSA_ana', '-f', app.FILES.prefix + 'info'])
        if g.wait():
            end = 1
    if end:
        logging.error('Unable to start gmx_MMPBSA_ana...')
    logging.info('Finalized...')
    return end

def gmxmmpbsa(argString: str):
    argv = shlex.split(argString)
    logging.basicConfig(
        level=logging.INFO,
        format="[%(levelname)-7s] %(message)s",
        handlers=[
            logging.FileHandler("gmx_MMPBSA.log", 'w'),
            logging.StreamHandler()])
    # Just for compatibility as mpi4py works as serial when run without mpirun
    # (since v1.4.2)
    if len(argv) > 1 and argv[1] in ['MPI', 'mpi']:
        args = argv
        args.pop(1)  # remove mpi arg before passed it to app
        from mpi4py import MPI
        # try:
        # except ImportError:
        #     GMXMMPBSA_ERROR('Could not import mpi4py package! Use serial version or install mpi4py.')
    else:
        # If we're not running "gmx_MMPBSA MPI", bring MPI into the top-level namespace
        # (which will overwrite the MPI from mpi4py, which we *want* to do in serial)
        from GMXMMPBSA.fake_mpi import MPI
        args = argv
    # Set up error/signal handlers
    pbsaMain.setup_run()

    # Instantiate the main MMPBSA_App
    app = pbsaMain.MMPBSA_App(MPI)

    # Read the command-line arguments
    try:
        app.get_cl_args(args[1:])
    except CommandlineError as e:
        sys.stderr.write('%s: %s' % (type(e).__name__, e) + '\n')
        sys.exit(1)

    # Perform our MMPBSA --clean now
    if app.FILES.clean:
        sys.stdout.write('Cleaning temporary files and quitting.\n')
        app.remove(0)
        return 0

    # See if we wanted to print out our input file options
    if app.FILES.infilehelp:
        app.input_file.print_contents(sys.stdout)
        sys.exit(0)

    # If we're not rewriting output do whole shebang, otherwise load info and parms
    # Throw up a barrier before and after running the actual calcs
    if not app.FILES.rewrite_output:
        try:
            app.read_input_file()
        except InputError as e:
            sys.stderr.write('%s: %s' % (type(e).__name__, e) + '\n')
            sys.stderr.write('  Enter `%s --help` for help\n' %
                             (split(sys.argv[0])[1]))
            sys.exit(1)
        app.process_input()
        app.check_for_bad_input()
        app.make_prmtops()
        app.loadcheck_prmtops()
        app.file_setup()
        app.run_mmpbsa()
    # If we are rewriting output, load the info and check prmtops
    else:
        info = InfoFile(app)
        info.read_info()
        app.make_prmtops()
        app.loadcheck_prmtops()

    # Now we parse the output, print, and finish
    app.parse_output_files()
    app.write_final_outputs()
    #app.finalize()
    exitCode = finalize(app)
    if exitCode != 0:
        raise Exception('ERROR finalize the run with exit code %d'%exitCode)

def obtain_number_of_traj_frame(topfile, trajfile):
    u = mda.Universe(topfile)
    u.load_new(trajfile)
    return len(u.trajectory)

def mmpbsa_pdb(pdbfile, groupIDs, outdir, prep=True, DEBUG=False):
    pdbfile = os.path.abspath(pdbfile)
    pwd = os.getcwd()
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    rundir = os.path.join(outdir, 'mmpbsa')
    if not os.path.exists(rundir):
        os.mkdir(rundir)
    os.chdir(rundir)
    if prep:
        outfile = os.path.split(pdbfile)[-1]#[:-4] + "_prep.pdb"
        pdb2pqr(pdbfile, outfile)
        pdbfile = outfile
    indexFile = 'index.ndx'
    gmx_index_group(pdbfile, selectedGroup=groupIDs, indexFile=indexFile)
    pbsaFile = 'mmpbsa.in'
    generate_input_file(outfile=pbsaFile)
    select_group(pdbfile, [groupIDs[1]], 'ligand.pdb')
    cmd = 'antechamber -i {input} -fi pdb -o {output} -fo mol2 -c gas'.format(input='ligand.pdb', output='ligand.mol2')
    RC = os.system(cmd)
    if RC != 0:
        raise Exception('ERROR run command %s'%cmd)
    receptorIndex, ligandIndex = 1,2
    parasDict = {
        "pbsaFile" : pbsaFile,
        "topFile"  : pdbfile,
        'indexFile': indexFile,
        'receptorIndex': receptorIndex,
        'ligandIndex': ligandIndex,
        'trajFile' : pdbfile,
        'ligandmol2':'ligand.mol2',
    }
    argString = '_ -O -i {pbsaFile} -cs {topFile} -ci {indexFile} -cg {receptorIndex} {ligandIndex} -ct {trajFile} -nogui -lm {ligandmol2}'.format(**parasDict)
    print(argString)
    gmxmmpbsa(argString)
    outfile = 'FINAL_RESULTS_MMPBSA.dat'
    if not os.path.exists(outfile):
        raise Exception('ERROR: Not found the result file %s'%outfile)
    cmd = 'cp %s ../'%(outfile)
    os.system(cmd)
    os.chdir(pwd)
    if not DEBUG:
        cmd = 'rm -rf %s'% rundir
        print('Clean output')
        RC = os.system(cmd)
        if RC != 0:
            raise Exception('ERROR run command %s'%cmd)

def mmbpsa_traj(tprfile, trajfile, groupIDs, outdir, indexFile=None, DEBUG=False):
    '''
    1. generate an index file for the selected group
        2. select the group from the tpr file
        3. generate a mol2 file for the ligand
        4. run mmpbsa
        5. copy the result file to the output directory
        6. clean the output directory
    
    # Python
    def mmpbsa_traj_with_index(tprfile, trajfile, groupIDs, outdir, indexFile=None, DEBUG=False):
        tprfile = os.path.abspath(tprfile)
        trajfile = os.path.abspath(trajfile)
        nFrame = obtain_number_of_traj_frame(tprfile, trajfile)
        pwd = os.getcwd()
        if not os.path.exists(outdir):
            os.mkdir(
    
    Args:
      tprfile: The input tpr file
      trajfile: the trajectory file
      groupIDs: a list of two group IDs, e.g. [1,2]
      outdir: the directory where the output files will be written
      indexFile: The index file that contains the group definitions.
      DEBUG: If you want to debug the script, set this to True. Defaults to False
    '''
    tprfile = os.path.abspath(tprfile)
    trajfile = os.path.abspath(trajfile)
    nFrame = obtain_number_of_traj_frame(tprfile, trajfile)
    pwd = os.getcwd()
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    rundir = os.path.join(outdir, 'mmpbsa')
    if not os.path.exists(rundir):
        os.mkdir(rundir)
    os.chdir(rundir)
    if indexFile is None:
        indexFile = 'index.ndx'
        gmx_index_group(tprfile, selectedGroup=groupIDs, indexFile=indexFile)
    pbsaFile = 'mmpbsa.in'
    generate_input_file(outfile=pbsaFile, endFrame=nFrame)
    select_group(tprfile, [groupIDs[1]], 'ligand.pdb', trajfile)
    cmd = 'antechamber -i {input} -fi pdb -o {output} -fo mol2 -c gas'.format(input='ligand.pdb', output='ligand.mol2')
    RC = os.system(cmd)
    if RC != 0:
        raise Exception('ERROR run command %s'%cmd)
    receptorIndex, ligandIndex = 1,2
    parasDict = {
        "pbsaFile" : pbsaFile,
        "topFile"  : tprfile,
        'indexFile': indexFile,
        'receptorIndex': receptorIndex,
        'ligandIndex': ligandIndex,
        'trajFile' : trajfile,
        'ligandmol2':'ligand.mol2',
    }
    argString = '_ -O -i {pbsaFile} -cs {topFile} -ci {indexFile} -cg {receptorIndex} {ligandIndex} -ct {trajFile} -nogui -lm {ligandmol2}'.format(**parasDict)
    print(argString)
    gmxmmpbsa(argString)
    outfile = 'FINAL_RESULTS_MMPBSA.dat'
    if not os.path.exists(outfile):
        raise Exception('ERROR: Not found the result file %s'%outfile)
    cmd = 'cp %s ../'%(outfile)
    os.system(cmd)
    os.chdir(pwd)
    if not DEBUG:
        cmd = 'rm -rf %s'% rundir
        print('Clean output')
        RC = os.system(cmd)
        if RC != 0:
            raise Exception('ERROR run command %s'%cmd)

def main():
    parser = argparse.ArgumentParser(description='Free energy calcaulated by MMPBSA method.')
    parser.add_argument('-i', dest='INP', help='A pdb file or a tpr file to calculate the free energy.', required=True)
    parser.add_argument('-t', dest='TRAJ', help='A trajectory file contains many structure frames. File format: xtc, pdb, gro...', default=None)
    parser.add_argument('-o', dest='OUTP', help='Output floder to save results.', default=None)
    parser.add_argument('-ndx', dest='ndx', help='Index file.', default=None)
    parser.add_argument('-D', dest='DEBUG', help='DEBUG model, keep all the files.', default=False, action='store_true')

    args = parser.parse_args()
    print(args)
    #exit(0)
    pdbfile, trajfile, outdir, debug, indexFile = args.INP, args.TRAJ, args.OUTP, args.DEBUG, args.ndx
    doc = 'ERROR: trajfile reuired for input a tpr file.'
    ext = pdbfile.split('.')[-1]
    if ext=='tpr' and trajfile is None:
        print(doc)
    ## 1 Select Groups
    groups = detect_group(pdbfile)
    print('='*80)
    print("Group list: %s"%pdbfile)
    print_group(groups)
    inputstr = input('\nPlease select to group as receptor and ligand(eg 1,2):')
    try:
        groupIDs = np.array(list(map(int, inputstr.split(',')))) - 1
    except:
        raise Exception('ERROR: You must input as 1,2!')
    ## 
    if trajfile:
        mmbpsa_traj(pdbfile, trajfile, groupIDs, outdir, indexFile, DEBUG=debug)
    else:
        #groupIDs = [24,2]
        mmpbsa_pdb(pdbfile, groupIDs, outdir, DEBUG=debug)

if __name__ == "__main__":
    main()
    
