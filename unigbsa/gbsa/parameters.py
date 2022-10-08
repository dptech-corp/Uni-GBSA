from unigbsa.settings import PBSA_VERSION, PBSA_PARAMETER_FILE
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

def generate_input_file_v143(pbsaParas, outfile='mmpbsa.in') -> None:
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
    mode = pbsaParas['mode'] 
    name = ''
    modes = ['gb', 'pb', 'gb+pb', 'pb+gb']
    if mode not in modes:
        raise Exception('Unknown type of mode: %s'%mode)
    args = {
        'sysName':name+"_"+mode,
        'mode':mode,
        "startFrame": pbsaParas['startFrame'] if 'startFrame' in pbsaParas  else 1,
        "endFrame": pbsaParas['endFrame'] if 'endFrame' in pbsaParas  else 1,
        "interval": pbsaParas['interval'] if 'interval' in pbsaParas  else 1,
        "temperature": pbsaParas['temperature'] if 'temperature' in pbsaParas  else 300,
        "igbValue": pbsaParas['igbValue'] if 'igbValue' in pbsaParas  else 2,
        "indi": pbsaParas['indi'] if 'indi' in pbsaParas  else 1.0,
        "exdi": pbsaParas['exdi'] if 'exdi' in pbsaParas  else 80.0,
    }
    line = generalstring.format(**args)
    if 'gb' in mode:
        line += gbstring.format(**args)
    if 'pb' in mode:
        line += pbstring.format(**args)
    if "decompose" in pbsaParas and pbsaParas['decompose']:
        line += decostring
    with open(outfile, 'w') as fw:
        fw.write(line)
    return outfile

def generate_input_file_v152(pbsaParas, outfile='mmpbsa.in'):
    modes = pbsaParas['modes'].split(',') + ['general']
    if 'indi' in pbsaParas:
        pbsaParas["intdiel"] = pbsaParas['indi']
    if 'exdi' in pbsaParas:
        pbsaParas['extdiel'] = pbsaParas['exdi']
    with open (PBSA_PARAMETER_FILE) as fr:
        lines = fr.readlines()
    mode = None
    MF = False
    with open(outfile, 'w') as fw:
        for line in lines:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            if line.startswith('&'):
                if MF:
                    fw.write('/\n\n')
                    MF = False
                mode = line.split('&')[1]
                if mode in modes:
                    fw.write(line+"\n")
                    MF = True
            elif mode in modes and '=' in line:
                lineList = line.split('#')
                comment = lineList[-1]
                varlist = lineList[0].split('=')
                varname = varlist[0].strip()
                varvalue = varlist[1].strip()
                if varname in pbsaParas:
                    if '"' in varvalue or "'" in varvalue:
                        varvalue = '"%s"'%pbsaParas[varname].replace('"','')
                    else:
                        varvalue = pbsaParas[varname]
                line = '  %-21s= %-47s# %s\n'%(varname, varvalue, comment)
                fw.write(line)
        if MF:
            fw.write('/\n')
    return outfile

def set_parameters(mmpbsafile, key, value):
    with open(mmpbsafile) as fr:
        lines = fr.readlines()
    with open(mmpbsafile, 'w') as fw:
        for line in lines:
            if key in line:
                lineList = line.strip().split(key)
                if ',' in lineList[1]:
                    tmp = lineList[1].split(',')
                    line = lineList[0] + key +'=%s'%str(value) + ','.join(tmp[1:]) + '\n'
                else:
                    line = lineList[0] + key + '=%s'%str(value) + '\n'
            fw.write(line)

if PBSA_VERSION >=1.5:
    generate_input_file = generate_input_file_v152
else:
    generate_input_file = generate_input_file_v143

if __name__ == "__main__":
    pass
    #generate_input_file()

