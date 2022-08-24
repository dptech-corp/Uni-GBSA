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
#        'mdanalysis>=2.0.0',
#        'pdb2pqr',
]
setup(
    name = 'hmtpbsa',
    version=get_version("hmtpbsa/version.py"),
    author='dptech.net',
    author_email='hermite@dptech.net',
    description=('MMBPSA tools for calculate energy.'),
    url='https://github.com/dptech-corp/Hermite-MMPBSA',
    license=None,
    keywords='MMPBSA',
 #   install_requires=install_requires,
    packages=find_packages(),
    zip_safe = False,
    #packages=packages,
    entry_points={'console_scripts': [
         'hmtpbsa-pipeline = hmtpbsa.pipeline:main',
         'hmtpbsa-traj = hmtpbsa.pbsa.pbsarun:main',
         'hmtpbsa-pbc = hmtpbsa.CLI:PBC_main',
         'hmtpbsa-buildtop = hmtpbsa.CLI:topol_builder',
         'hmtpbsa-buildsys = hmtpbsa.CLI:simulation_builder',
         'hmtpbsa-md = hmtpbsa.CLI:simulation_run',
     ]},
    include_package_data=True
)
