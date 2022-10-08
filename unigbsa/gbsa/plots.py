import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl

from unigbsa.gbsa import io

fontdir = Path(__file__).absolute().parent.parent.joinpath('data/Calibri')
fontfiles = mpl.font_manager.findSystemFonts(str(fontdir))
for fontfile in fontfiles:
    mpl.font_manager.fontManager.addfont(fontfile)
fontdict = {
    'font.family' : 'Calibri',
    'font.size' : 12,
    'figure.max_open_warning' : 1000
}

DPI=150
mpl.rcParams.update(fontdict)

def plot_deltaG_component(df, type='GB', ax=None, outfile=None, dpi=DPI):
    color = {
        'GB': ['green']*7 + ['blue']*2 + ['green', 'blue', 'red'],
        'PB': ['green']*7 + ['blue']*3 + ['green', 'blue', 'red']
    }
    mask = (df['type']==type) & (df['group']=='Delta (Complex - Receptor - Ligand)')
    data = df[mask].reindex()
    x = data['component'].apply(lambda x: str(x).replace('ΔTOTAL', 'ΔG')).values
    x = [r'$\mathregular{\delta G}$']
    
    if ax is None:
        fig, ax = plt.subplots()
    ax.bar(x=data['component'], height=data['average'], yerr=data['SD'], color=color[type])
    xmin, xmax = ax.get_xlim()
    ax.hlines(y=0, xmin=xmin, xmax=xmax, color='black', linewidth=1)
    ax.set_ylabel('energy (kcal/mol)')
    for xtick in ax.get_xticklabels():
        xtick.set_rotation(50)
    ax.set_xlim(xmin, xmax)
    fig.tight_layout()
    if outfile:
        fig.savefig(outfile, dpi=dpi)
    plt.close(fig)

def plot_deltaG_DECOMP(df, chain=None, key='TOTAL', outfile=None, dpi=DPI): 
    df['restype'] = df['residue'].apply(lambda x: str(x)[0])
    df['resname'] = df['residue'].apply(lambda x: str(x)[3:])
    data = df[df['restype']=='R']
    y = data[key].values.astype('float')
    yerr = data[key+' Std. Err. of Mean'].values.astype('float')

    fig, ax = plt.subplots()
    ax.bar(x=data['resname'], height=y, yerr=yerr, color='blue')
    ax.set_ylabel('energy (kcal/mol)')
    for xtick in ax.get_xticklabels():
        xtick.set_rotation(50)
    if key == 'TOTAL':
        key = r'$\mathregular {\Delta G}$'
    ax.set_title(key)
    fig.tight_layout()
    if outfile:
        fig.savefig(outfile, dpi=dpi)
    plt.close(fig)

def plot_res_deltaG_component(df, resname, chain=None, outfile=None, dpi=DPI): 
    #df['restype'] = df['residue'].apply(lambda x: str(x)[0])
    #df['resname'] = df['residue'].apply(lambda x: str(x)[3:])
    mask = (df['restype']=='R') & (df['residue']==resname)
    data = df[mask]
    if len(data)==0:
        return
    keys = ['Internal','van der Waals','Electrostatic','Polar Solvation', 'Non-Polar Solv.', 'TOTAL']
    keys_std = [ k+' Std. Err. of Mean' for k in keys ]
    labels = ['Internal','van der Waals','Electrostatic','Polar Solv.', 'Non-Polar Solv.', 'TOTAL']
    y = data[keys].values.astype('float')[0]
    yerr = data[keys_std].values.astype('float')[0]
    fig, ax = plt.subplots()
    ax.bar(x=keys, height=y, yerr=yerr, color='blue')
    ax.set_xticks(keys, labels)
    ax.set_ylabel('energy (kcal/mol)')
    for xtick in ax.get_xticklabels():
        xtick.set_rotation(50)
    ax.set_title(resname)
    fig.tight_layout()
    if outfile:
        fig.savefig(outfile, dpi=dpi)
    plt.close(fig)

def plot_deltaG_traj(dic, key='TOTAL', lw=3, outfile=None, dpi=DPI):
    color = {
        'GB':'red',
        'PB': 'blue'
    }
    title = r'$\mathregular{\Delta G_{%s}}$'%key[1:] if key.startswith('G') else r'$\mathregular{\Delta G_{%s}}$'%key
    fig, ax = plt.subplots()
    for ptype in dic.keys():
        df = dic[ptype]
        x = df['Frame #'].values.astype(int)
        if key in df:
            y = df[key].values.astype(float)
            ax.plot(x-1, y, color=color[ptype], linewidth=lw, label=ptype)
    ax.set_xlabel('# frame')
    ax.set_ylabel('energy (kcal/mol)')
    ax.set_title(title)
    fig.legend(loc=(0.15,0.75))
    fig.tight_layout()
    if outfile:
        fig.savefig(outfile, dpi=dpi)
    plt.close(fig)

def plot_deltaG_traj_decomp(dic, res, key='TOTAL', lw=3, outfile=None, dpi=DPI):
    color = {
        'GB':'red',
        'PB': 'blue'
    }
    fig, ax = plt.subplots()
    for ptype in dic.keys():
        mask = dic[ptype]['Residue'] == res
        df = dic[ptype][mask]
        x = df['Frame #'].values.astype(int)
        if key in df:
            y = df[key].values.astype(float)
            ax.plot(x-1, y, color=color[ptype], linewidth=lw, label=ptype)
    ax.set_xlabel('# frame')
    ax.set_ylabel('energy (kcal/mol)')
    title = res[3:] + ' ' + key
    ax.set_title(title)
    fig.legend(loc=(0.15,0.75))
    fig.tight_layout()
    if outfile:
        fig.savefig(outfile, dpi=dpi)
    plt.close(fig)
    

def analysis_FINAL(finalfile='FINAL_RESULTS_MMPBSA.dat', outdir='analysis'):
    df = io.read_FINAL_output(finalfile)
    cwd = os.getcwd()
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    os.chdir(outdir)
    plot_deltaG_component(df, type='GB', outfile='GB-deltaG.png')
    plot_deltaG_component(df, type='PB', outfile='PB-deltaG.png')
    df.to_csv('deltaG.csv', index=False)
    os.chdir(cwd)

def analysis_DECOMP(decompfile='FINAL_DECOMP_MMPBSA.dat', outdir='analysis'):
    dic = io.read_DECOMP_output(decompfile)
    cwd = os.getcwd()
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    os.chdir(outdir)
    decomdir = 'decomp'
    if not os.path.exists(decomdir):
        os.mkdir(decomdir)
    keys = list(dic.keys())
    residues = dic[keys[0]]['residue']
    for key in keys:
        plot_deltaG_DECOMP(dic[key], outfile='%s-deltaG-decomp.png'%key)
        for residue in residues:
            outfile = os.path.join(decomdir, '%s-deltaG-decomp%s.png'%(key, residue.replace(':','_')))
            plot_res_deltaG_component(dic[key], residue, outfile=outfile)
        dic[key].to_csv('%s-deltaG-decomp.csv'%key, index=False)
    os.chdir(cwd)

def analysis_traj_EO(dofile='EO.csv', outdir='analysis'):
    dic = io.read_EO_output(dofile)
    cwd = os.getcwd()
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    os.chdir(outdir)
    trajdir = 'traj'
    if not os.path.exists(trajdir):
        os.mkdir(trajdir)
    s = 'BOND\tANGLE\tDIHED\tVDWAALS\tEEL\t1-4 VDW\t1-4 EEL\tEGB\tESURF\tGGAS\tGSOLV\tTOTAL'
    keys = s.split('\t')

    for k in keys:
        outfile = os.path.join(trajdir, 'traj_decom_%s.png'%k.replace(' ', '_'))
        if k=='TOTAL':
            outfile = 'traj_deltaG.png'
        plot_deltaG_traj(dic, key=k, outfile=outfile)
    for key in dic.keys():
        dic[key].to_csv('EO_%s.csv'%key, index=False)
    os.chdir(cwd)

def analysis_traj_DEO(deofile='DEO.csv', outdir='analysis'):
    dic = io.read_DEO_output(deofile)
    cwd = os.getcwd()
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    os.chdir(outdir)
    trajdir = 'traj'
    if not os.path.exists(trajdir):
        os.mkdir(trajdir)
    keys = list(dic.keys())

    residues = dic[keys[0]]['Residue'].unique()
    for residue in residues:
        outfile = os.path.join(trajdir, 'traj_res_%s.png'%(residue.replace(':', '_')))
        plot_deltaG_traj_decomp(dic, residue, outfile=outfile)
    for key in keys:
        dic[key].to_csv('DEO_%s.csv'%key, index=False)
    os.chdir(cwd)
