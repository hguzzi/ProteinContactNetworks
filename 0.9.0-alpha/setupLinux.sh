#!/bin/bash
echo Setup Conda Environment
if [[ $OSTYPE == 'linux-gnu'* ]]; then
echo 'Linux'
	conda create -n PCN python=3.8.3
	conda activate PCN 
	sudo apt install pip
	sudo apt install git
	pip install scipy
	pip install ProDy
	pip install numpy
	pip install fuzzy-c-means
	pip install git+https://github.com/palash1992/GEM.git
	pip install git+https://github.com/GiulioRossetti/cdlib.git
	conda install -c schrodinger pymol
	pip install python-dateutil
	pip install pytz
	python create_config_file.py
fi
