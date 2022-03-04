'''
see the https://valdes-tresanco-ms.github.io/gmx_MMPBSA/input_file/#the-input-file
'''

from quopri import decodestring


generalstring = '''

&general
sys_name = {sysName}, startframe={startFrame}, endframe={endFrame}, verbose=2, interval={interval}, temperature={temperature},
protein_forcefield="oldff/leaprc.ff99SB",
/
'''

gbstring = '''
&gb
igb={igbValue}, saltcon=0.150
/
'''

pbstring = '''
&pb
istrng=0.15, fillratio=4.0
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

def generate_input_file(mode='gb', outfile='mmpbsa.in', startFrame=1, endFrame=1, interval=1, temperature=300, igbValue=2, name='Calculate', decompose=False) -> None:

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
    }
    line = generalstring.format(**args)
    if 'gb' in mode:
        line += gbstring.format(**args)
    if 'pb' in mode:
        line += pbstring
    if decompose:
        line += decostring
    with open(outfile, 'w') as fw:
        fw.write(line)
    return outfile


if __name__ == "__main__":
    pass
    #generate_input_file()

