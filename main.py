from GMXMMPBSA.commandlineparser import index, topology
from groups import gmx_index_group
import sys
from os.path import split
import logging
import shlex


from parameters import generate_input_file

try:
    from GMXMMPBSA.exceptions import InputError, CommandlineError, GMXMMPBSA_WARNING
    from GMXMMPBSA.infofile import InfoFile
    from GMXMMPBSA import main
    from GMXMMPBSA.commandlineparser import anaparser, testparser
except ImportError:
    import os
    amberhome = os.getenv('AMBERHOME') or '$AMBERHOME'
    raise ImportError('Could not import Amber Python modules. Please make sure '
                      'you have sourced %s/amber.sh (if you are using sh/ksh/'
                      'bash/zsh) or %s/amber.csh (if you are using csh/tcsh)' %
                      (amberhome, amberhome))
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
    main.setup_run()

    # Instantiate the main MMPBSA_App
    app = main.MMPBSA_App(MPI)

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
        sys.exit(0)

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
    app.finalize()

def mmpbsa_pdb():
    pdbfile = 'test_prep.pdb'
    indexFile = 'index.ndx'
    gmx_index_group(pdbfile, indexFile=indexFile)
    pbsaFile = 'mmpbsa.in'
    generate_input_file(outfile=pbsaFile)
    receptorIndex, ligandIndex = 1,2
    parasDict = {
        "pbsaFile" : pbsaFile,
        "topFile"  : pdbfile,
        'indexFile': indexFile,
        'receptorIndex': receptorIndex,
        'ligandIndex': ligandIndex,
        'trajFile' : pdbfile
    }
    argString = '_ -O -i {pbsaFile} -cs {topFile} -ci {indexFile} -cg {receptorIndex} {ligandIndex} -ct {trajFile} -nogui'.format(**parasDict)

    gmxmmpbsa(argString)

def mmbpsa_traj():
    tprfile = 'examples/Protein_ligand/ST/com.tpr'
    trajfile = 'examples/Protein_ligand/ST/com_traj.xtc'



def test_group_select():
    from groups import detect_group, select_group, print_group, gmx_index_group
    from prepare import pdb2pqr
    pdbfile = './example/2fvy.pdb'
    pdbfile = '../example/1ycr.pdb'
    pdbfile = '../com_traj.pdb'
    
    groups = detect_group(pdbfile)
    print_group(groups)
    select_group(pdbfile, [3,4], outfile='test.pdb')
    
    outfile = 'test_prep.pdb'
    pdb2pqr('test.pdb', outfile)
    gmx_index_group(outfile)


def gen_index():
    pdbfile = 'test_out.pdb'
    gmx_index_group(pdbfile)

if __name__ == "__main__":
    #print(sys.argv)
    #gmxmmpbsa()
    test_group_select()
    #gen_index()
    #mmpbsa_pdb()
