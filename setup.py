import codecs
import os.path
from setuptools import setup, find_packages

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


install_requires = [
        'acpype',
        'gmx_MMPBSA>=1.5.6',
        'lickit'
]
setup(
    name = 'unigbsa',
    version=get_version("unigbsa/version.py"),
    author='dptech.net',
    author_email='hermite@dptech.net',
    description=('MMPB(GB)SA tools for calculate energy.'),
    url='https://github.com/dptech-corp/Uni-GBSA',
    license=None,
    keywords='MMPBSA MMGBSA',
    install_requires=install_requires,
    packages=find_packages(),
    zip_safe = False,
    #packages=packages,
    entry_points={'console_scripts': [
         'unigbsa-pipeline = unigbsa.pipeline:main',
         'unigbsa-traj = unigbsa.CLI:traj_pipeline',
         'unigbsa-pbc = unigbsa.CLI:PBC_main',
         'unigbsa-buildtop = unigbsa.CLI:topol_builder',
         'unigbsa-buildsys = unigbsa.CLI:simulation_builder',
         'unigbsa-md = unigbsa.CLI:simulation_run',
         'unigbsa-plots = unigbsa.CLI:mmpbsa_plot'
     ]},
    include_package_data=True
)
