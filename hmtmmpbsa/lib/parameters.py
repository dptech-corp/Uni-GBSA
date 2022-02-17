'''
see the https://valdes-tresanco-ms.github.io/gmx_MMPBSA/input_file/#the-input-file
'''

stability = '''
Sample input file for GB calculation
This input file is meant to show only that gmx_MMPBSA works. Althought,
we tried to used the input files as recommended in the Amber manual,
some parameters have been changed to perform more expensive calculations
in a reasonable amount of time. Feel free to change the parameters 
according to what is better for your system.

&general
startframe={startFrame}, endframe={endFrame}, verbose=2, interval={interval}, temperature={temperature}
protein_forcefield="oldff/leaprc.ff99SB",
/

&gb
igb={igbValue}, saltcon=0.150,
/
'''



def generate_input_file(mode=1, outfile='mmpbsa.in', startFrame=1, endFrame=1, interval=1, temperature=300, model='gb') -> None:
    """[summary]

    Args:
        mode (int, optional): Calculate mode: 
                   1: protein-protein, 
                   2: protein-ligand. 
                     Defaults to 1.
        outfile (str, optional): [description]. Defaults to 'mmpbsa.in'.
        startFrame (int, optional): [description]. Defaults to 1.
        endFrame (int, optional): [description]. Defaults to 1.
        interval (int, optional): [description]. Defaults to 1.
        temperature (int, optional): [description]. Defaults to 300.
    """
    igbValue = {
        1: 2, # protein-protein
        2: 5  # protein-ligand
    }
    if mode not in igbValue.keys():
        raise Exception('suport value for mode: %d'%mode)
    args = {
        "startFrame":startFrame,
        "endFrame": endFrame,
        "interval": interval,
        "temperature": temperature,
        "igbValue": igbValue[mode]
    }
    line = stability.format(**args)
    with open(outfile, 'w') as fw:
        fw.write(line)

if __name__ == "__main__":
    pass
    #generate_input_file()

