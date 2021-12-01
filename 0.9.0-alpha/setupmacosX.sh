#!/bin/bash
echo Setup Conda Environment
if [[ $OSTYPE == 'darwin'* ]]; then
  echo 'macOS'
  conda create -n PCN python=3.8.3
  conda activate PCN 
  pip install scipy
  pip install ProDy
  pip install numpy
  pip install scipy
  pip install fuzzy-c-means
  pip install git+https://github.com/palash1992/GEM.git
  pip install git+https://github.com/GiulioRossetti/cdlib.git
  pip install python-dateutil
  pip install pytz
  conda install -c schrodinger pymol
fi




