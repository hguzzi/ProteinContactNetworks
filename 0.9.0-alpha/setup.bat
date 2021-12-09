call activate.bat 
call conda create -n PCN python=3.8.3
call conda install -n PCN pip
call python create_config_file.py
pause