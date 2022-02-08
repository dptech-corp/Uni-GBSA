from setuptools import setup, find_packages

install_requires = [
#        'mdanalysis>=2.0.0',
#        'pdb2pqr',
]
setup(
    name = 'hmtmmpbsa',
    version='0.0.1',
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
         'hmtmmpbsa = hmtmmpbsa.main:main',
     ]},
    include_package_data=True
)
