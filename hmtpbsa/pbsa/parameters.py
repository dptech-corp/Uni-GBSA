'''
see the https://valdes-tresanco-ms.github.io/gmx_MMPBSA/input_file/#the-input-file
'''

generalstring = '''

&general
sys_name = {sysName}, startframe={startFrame}, endframe={endFrame}, verbose=2, interval={interval}, temperature={temperature}
protein_forcefield="oldff/leaprc.ff99SB",
/
'''

gbstring = '''
&gb
igb={igbValue}, saltcon=0.150, intdiel={indi}
/
'''

pbstring = '''
&pb
istrng=0.15, fillratio=4.0, indi={indi}, exdi={exdi}
/
'''

decostring = '''
&decomp
idecomp=2, dec_verbose=0,
# This will print all residues that are less than 4 Å between
# the receptor and the ligand
print_res="within 5"
/
'''

def generate_input_file(mode='gb', outfile='mmpbsa.in', startFrame=1, endFrame=1, interval=1, temperature=300, igbValue=2, name='Calculate', decompose=False, indi=1.0, exdi=80.0) -> None:
    """
    Generate a mmpbsa.in file for a given mode
    
    Args:
      mode: gb, pb, gb+pb, pb+gb. Defaults to gb
      outfile: the name of the input file. Defaults to mmpbsa.in
      startFrame: the first frame to be used in the calculation. Defaults to 1
      endFrame: The last frame of the trajectory to be used. Defaults to 1
      interval: the interval between frames to be used in the calculation. Defaults to 1
      temperature: 300. Defaults to 300
      igbValue: GB model. 1:igb=5; 2:igb=7; 3:igb=8. Defaults to 2
      name: the name of the system. Defaults to Calculate
      decompose: whether to decompose the GB/PB energy into individual contributions. Defaults to False
      indi: the initial distance between the two ligand atoms in the GB calculation
      exdi: The end of the dielectric region.
    
    Returns:
      The name of the input file.
    """
    modes = ['gb', 'pb', 'gb+pb', 'pb+gb']
    if mode not in modes:
        raise Exception('Unknown type of mode: %s'%mode)
    args = {
        'sysName':name+"_"+mode,
        'mode':mode,
        "startFrame":startFrame,
        "endFrame": endFrame,
        "interval": interval,
        "temperature": temperature,
        "igbValue": igbValue,
        "indi": indi,
        "exdi": exdi,
    }
    line = generalstring.format(**args)
    if 'gb' in mode:
        line += gbstring.format(**args)
    if 'pb' in mode:
        line += pbstring.format(**args)
    if decompose:
        line += decostring
    with open(outfile, 'w') as fw:
        fw.write(line)
    return outfile


if __name__ == "__main__":
    pass
    #generate_input_file()

