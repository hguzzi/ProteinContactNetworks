#!/bin/bash

if [[ $OSTYPE == 'linux-gnu'* ]] || [[ $OSTYPE == 'darwin'* ]];
then  

	echo $OSTYPE	
	
	echo "Preparing a new conda enviroment called PCN"
	conda create -n PCN python=3.8.3
	conda activate PCN 
		
	python -m ensurepip --upgrade
	pip install configparser
	cd pcn
	cd tools
	python create_config_file.py
	
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
	pip install pyqt5
fi

