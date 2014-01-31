from setuptools import setup
from glob import glob

setup(name='NGTS_workpackage',
        version='0.0.1',
        description='NGTS photometry workpackages',
        author='Tom Louden',
        author_email='t.m.louden@warwick.ac.uk',
        maintainer='Simon Walker',
        maintainer_email='s.r.walker101@googlemail.com',
        url='http://github.com/tomlouden/NGTS_ZLP',
        packages=['NGTS_workpackage', ],
        scripts=glob('bin/*.py'),
        long_description=open('README.md').read(),
        install_requires=['astropy>=0.3',
            'docopt',
            'numpy>=1.8',
            'scipy',
            ]
        )
