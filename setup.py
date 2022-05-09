from setuptools import setup, find_packages

install_requires = [
#        'mdanalysis>=2.0.0',
#        'pdb2pqr',
]
setup(
    name = 'hmtpbsa',
    version='0.0.5_pre',
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
