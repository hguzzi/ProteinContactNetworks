@echo off

echo "Preparing a new conda enviroment called PCN"
call activate.bat 
call conda create -n PCN python=3.8.3
call conda install -n PCN pip
call pip install configparser
cd pcn
cd tools
call python create_config_file.py

echo "Installing libraries"
call conda activate PCN
call pip install numpy
call pip install regex
call pip install scipy
call pip install fuzzy-c-means
call pip install git+https://github.com/palash1992/GEM.git
call pip install git+https://github.com/GiulioRossetti/cdlib.git
call conda install -c schrodinger pymol
call pip install pytz
call pip install python-dateutil
call pip install git+https://github.com/eliorc/node2vec.git

pause