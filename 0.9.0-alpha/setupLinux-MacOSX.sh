#!/bin/bash

if [[ $OSTYPE == 'linux-gnu'* ]]||[[ $OSTYPE == 'darwin'* ]]; then
  
	echo $OSTYPE
	
	echo "Preparing a new conda enviroment called PCN"
	conda create -n PCN python=3.8.3
	conda activate PCN 
	sudo apt install pip
	sudo apt install git
	python create_config_file.py
	
	echo "Installing libraries"
	pip install scipy
	pip install numpy
	pip install fuzzy-c-means
	pip install git+https://github.com/palash1992/GEM.git
	pip install git+https://github.com/GiulioRossetti/cdlib.git
	conda install -c schrodinger pymol
	pip install python-dateutil
	pip install pytz
	python create_config_file.py

fi
