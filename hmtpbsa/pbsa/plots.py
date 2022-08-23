import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

fontfiles = mpl.font_manager.findSystemFonts('../data/Calibri/')
for fontfile in fontfiles:
    mpl.font_manager.fontManager.addfont(fontfile)
fontdict = {
    'font.family' : 'Calibri',
    'font.size' : 12
}

DPI=100
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

def plot_res_deltaG_component(df, resname, chain=None, outfile=None, dpi=DPI): 
    df['restype'] = df['residue'].apply(lambda x: str(x)[0])
    df['resname'] = df['residue'].apply(lambda x: str(x)[3:])
    mask = (df['restype']=='R') & (df['resname']==resname)
    data = df[mask]
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