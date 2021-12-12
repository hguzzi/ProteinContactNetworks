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
    
print('Config.ini file creation...')
print('Input the Directory in which you want to store the outputs')
print('The software will create three subdirectories')
output_path = str (input("Insert Root Path of the Outputs: "))
if(not os.path.isdir(output_path)):
    print("{} directory not exists.".format(output_path))
    end = False
    while(end): 
        output_choice = int(input("Press 0 if you want to create a new directory with path {} or type 1 if you want to retype the Output folder path")) 
        if(output_choice == 0):
            os.makedirs(output_path)
            end = True
        elif(output_choice == 1):
            output_path = str (input("Insert Root Path of the Outputs: "))     
            end = True
        #else: continue
            
print('Input the Directory containing Input Files')    
proteins_path = str (input("Insert Proteins filepath: "))
if(not os.path.isdir(proteins_path)):
    print("{} directory not exists.".format(proteins_path))
    end = False
    while(end): 
        proteins_choice = int(input("Press 0 if you want to create a new directory with path {} or type 1 if you want to retype the Protein folder path")) 
        if(proteins_choice == 0):
            os.makedirs(proteins_path)
            end = True
        elif(proteins_choice == 1):
            proteins_path = str (input("Insert Proteins filepath: "))      
        #else: continue
    
print('Please insert the path of the directory containing Adjacency Matrixs')
adj_filespath = str( input("Insert Adjacency matrix filepath: "))
if(not os.path.isdir(adj_filespath )):
    print("{} directory not exists.".format(adj_filespath ))
    end = False
    while(end): 
        adj_choice = int(input("Press 0 if you want to create a new directory with path {} or type 1 if you want to retype the adjacency matrixs folder path")) 
        if(adj_choice == 0):
            os.makedirs(adj_filespath)
            end = True
        elif(adj_choice == 1):
            adj_filespath = str( input("Insert Adjacency matrix filepath: "))   
        #else: continue
print('')


config = configparser.ConfigParser()
# Add the structure to the file we will create
config.add_section('user_paths')
config.set('user_paths', 'output_path', output_path)
config.set('user_paths', 'proteins_path', proteins_path)
config.set('user_paths', 'adj_filespath', adj_filespath)
# Write the new structure to the new file
with open(os.getcwd()+add_slash_to_path+"config.ini", 'w') as configfile:
    config.write(configfile)
