from setuptools import setup, find_packages

NAME = 'PCN-Miner'

DISTNAME = 'pcn'
MAINTAINER = 'Ugo Lomoio <ugo.lomoio@studenti.unicz.it>, Pietro Hiram Guzzi <hguzzi@unicz.it>'
DESCRIPTION = 'PCN-Miner: A tool for the analysis of Protein Contact Networks'
LONG_DESCRIPTION = open('README.md').read()
URL = 'https://github.com/hguzzi/ProteinContactNetworks'
KEYWORDS = ['Protein Contact Networks', 'Allostery', 'Proteins', 'Community Detection', 'Graph Embedding', 'Clustering', 'Spectral Clustering', 'Graphs']
LICENSE = 'CC0 1.0 Universal (CC0 1.0) Public Domain Dedication'
VERSION = '1.0.0'


INSTALL_REQUIRES = ['numpy','future','matplotlib','scikit-learn','tqdm','networkx>=2.4','demon','python-louvain>=0.16','nf1','scipy','pulp','seaborn','pandas',
                    'eva_lcd','bimlpa','markov_clustering','chinese_whispers','python-igraph','angel-cd','pooch','dynetx','thresholdclustering','pyclustering',
                    'cython','python-Levenshtein','regex','fuzzy-c-means','cdlib','pytz','python-dateutil','node2vec',
                    'GEM @ git+https://github.com/palash1992/GEM.git']

def setup_package():
    
    setup(
        name=NAME,
        version=VERSION,
        author=MAINTAINER,
        description=DESCRIPTION,
        url=URL,
        keywords=KEYWORDS,
        install_requires=INSTALL_REQUIRES,
        packages=find_packages(where=DISTNAME, exclude = 'gui'),
        package_dir = {'':'pcn'},
        license=LICENSE,
        long_description=LONG_DESCRIPTION,
        classifiers=['Intended Audience :: Science/Research',
                     'Intended Audience :: Developers',
                     'Development Status :: 4 - Beta',
                     'Operating System :: OS Independent',
                     'License :: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication',
                     'Programming Language :: Python :: 3',
                     ],
        )

if __name__ == "__main__":
    setup_package()