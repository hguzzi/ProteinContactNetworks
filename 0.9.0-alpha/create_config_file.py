import os 
from sys import platform
import configparser

if platform == "linux" or platform == "linux2":
    # linux
    add_slash_to_path = '/'
elif platform == "darwin":
    # OS X
    add_slash_to_path = '/'
elif platform == "win32":
    # Windows...
    add_slash_to_path = '\\' 
    
output_path = str (input("Insert Root Path of the Outputs: "))
proteins_path = str (input("Insert Proteins filepath: "))
adj_filespath = str( input("Insert Adjacency matrix filepath: "))


config = configparser.ConfigParser()
# Add the structure to the file we will create
config.add_section('user_paths')
config.set('user_paths', 'output_path', output_path)
config.set('user_paths', 'proteins_path', proteins_path)
config.set('user_paths', 'adj_filespath', adj_filespath)
# Write the new structure to the new file
with open(os.getcwd()+add_slash_to_path+"config.ini", 'w') as configfile:
    config.write(configfile)