from codecs import open
from os import path
from setuptools import setup, find_packages

here = path.abspath(path.dirname(__file__))
DISTNAME = 'src'
MAINTAINER = 'Ugo Lomoio'
DESCRIPTION = 'PCN-Miner: A tool for the analysis of Protein Contact Networks'
LONG_DESCRIPTION = open('README.md').read()
URL = 'https://github.com/hguzzi/ProteinContactNetworks'
KEYWORDS = ['Protein Contact Networks', 'Allostery', 'Proteins', 'Community Detection', 'Graph Embedding', 'Clustering', 'Spectral Clustering', 'Graphs']
LICENSE = 'CC0 1.0 Universal (CC0 1.0) Public Domain Dedication'
VERSION = '1.0.0'

with open(path.join(here, 'requirements.txt'), encoding='utf-8') as f:
    INSTALL_REQUIRES = f.read().splitlines()

def setup_package():
    
    setup(
        name=DISTNAME,
        version=VERSION,
        author=MAINTAINER,
        description=DESCRIPTION,
        url=URL,
        keywords=KEYWORDS,
        install_requires=INSTALL_REQUIRES,
        packages=find_packages(),
        package_dir={DISTNAME: 'pcn_miner'},
        license=LICENSE,
        long_description=LONG_DESCRIPTION,
        classifiers=['Intended Audience :: Science/Research',
                     'Intended Audience :: Developers',
                     'Development Status :: 4 - Beta',
                     'Operating System :: OS Independent',
                     'CC0 1.0 Universal (CC0 1.0) Public Domain Dedication',
                     'Programming Language :: Python :: 3',
                     ],
        )

if __name__ == "__main__":
    setup_package()
