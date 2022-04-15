#!/bin/bash


echo $OSTYPE	
	
		
echo "Installing libraries"
pip install scipy
pip install regex
pip install numpy
pip install fuzzy-c-means
pip install git+https://github.com/palash1992/GEM.git
pip install git+https://github.com/GiulioRossetti/cdlib.git
conda install -c schrodinger pymol 
pip install pytz
pip install python-dateutil
pip install git+https://github.com/eliorc/node2vec.git
pip install PyQt5
