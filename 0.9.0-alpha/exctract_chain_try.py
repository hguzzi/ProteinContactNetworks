import os
import networkx as nx
import pcn_final
import numpy as np
from scipy.linalg import eigh
import pcn_pymol_scripts
from sys import platform
import configparser
supported_algorithms_communities = ["louvain", "leiden", "walktrap", "asyn_fluidc", "greedy_modularity", "infomap", "spinglass"]


if platform == "linux" or platform == "linux2":
    # linux
    add_slash_to_path = '/'
    os.system('clear')
elif platform == "darwin":
    # OS X
    add_slash_to_path = '/'
    os.system('clear')
elif platform == "win32":
    # Windows...
    add_slash_to_path = '\\'
    supported_algorithms_communities.remove("infomap")

dict_supported_algorithms_communities = dict()    


dict_supported_algorithms_communities[str (0)] = "all"   
for i, algorithm in enumerate(supported_algorithms_communities):
    dict_supported_algorithms_communities[str (i+1)] = algorithm

   
if(os.path.isfile(os.getcwd()+add_slash_to_path+"config.ini")):
    print("config file found!")
    config_obj = configparser.ConfigParser()
    config_obj.read(os.getcwd()+add_slash_to_path+"config.ini")
    paths = config_obj["user_paths"]
    output_path = paths["output_path"]
    proteins_path = paths["proteins_path"]
    adj_filespath = paths["adj_filespath"]
    
    print("Paths in the config file: \n output path = {} \n proteins_path = {} \n adj_filespath = {}".format(output_path, proteins_path, adj_filespath))

else:#config file not found, create it
    print('No config file found... Start to create a new one')
    print('Input the Directory in which you want to store the outputs')
    print('The software will create three subdirectories')
    output_path = str (input("Insert Root Path of the Outputs: "))
    print('Input the Directory containing Input Files')    
    proteins_path = str (input("Insert Proteins filepath: "))
    print('Please insert the path of the directory containing Adjacency Matrixs')
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
        
is_dir_prot = os.path.isdir(proteins_path)
    
if is_dir_prot:
            
    if (not proteins_path.endswith(add_slash_to_path)):
        protein_list_dir = [file for file in os.listdir(proteins_path) if file.endswith('.pdb')]
        print("List of proteins in {}: {}".format(proteins_path, protein_list_dir))
            
        print('Please Insert Protein PDB Identifiers, separated by comma, without .pdb, e.g. 7nxc for 7nxc.pdb ')
        protein_choice = str(input("Enter proteins list items (splitted by ',') or 'all': "))
        
        if protein_choice.casefold() == 'all':
            proteins_list = protein_list_dir            
        else:
            proteins_list = [protein.casefold() for protein in protein_choice.replace(" ","").split(",")] #strip space
           
        proteins_path = proteins_path+add_slash_to_path
        pdb_list = proteins_list.copy()
        for i, protein in enumerate(pdb_list):
            if (not protein.endswith('.pdb')):
                pdb_list[i] = protein+'.pdb'
        pcn_final.checkIfFilesExists(pdb_list, "pdb", proteins_path) #add_slash_to_path is path separator             
    
    string_to_print = ""
    for i, algorithm in dict_supported_algorithms_communities.items():
        string_to_print =  string_to_print+"Type '{}' for '{}' algorithm/s \n".format(i, algorithm)
    
    algorithms_choice = []
    algorithms_choice_numeric = str( input("Choice one or a list (splitted by ',') Community Extraction algoritms from: \n{}Choice: ".format(string_to_print)))
    if(algorithms_choice_numeric == "0"):
        algorithms_choice = supported_algorithms_communities    
        
    if (algorithms_choice_numeric != "0"):
        if (isinstance(algorithms_choice_numeric, str)):
            if(algorithms_choice_numeric.split(',')):
                algorithms_choice_numeric =  [str(algorithm_choice_numeric) for algorithm_choice_numeric in algorithms_choice_numeric.replace(" ","").split(",") if(algorithm_choice_numeric != "0")]
        
        for algorithm_choice_numeric in algorithms_choice_numeric :
            algorithms_choice.append(dict_supported_algorithms_communities[algorithm_choice_numeric])            
        
        print('Selected algorithms: '+ str(algorithms_choice))
        
        not_supported_algorithms = []
        for algorithm_choice in algorithms_choice:
                      
            if (algorithm_choice not in supported_algorithms_communities):
                not_supported_algorithms.append(algorithm_choice)
                algorithms_choice.remove(algorithm_choice)
                        
        if (len(not_supported_algorithms)>0):
            print("Algorithms {} not supported yet.".format(str(not_supported_algorithms)))
            
            
    if ('asyn_fluidc' in algorithms_choice):
        k_choice = str(input("Entering k for Asyn FluidC: Enter an int, a list of ints (split with ','): "))
        if(k_choice.split(',')):
            ks =  [int(item) for item in k_choice.replace(" ","").split(",")]
            
min_ = float(input("Entering non covalent bonds threshold distance for PCN costruction: ") or 4.0)    
max_ = float(input("Entering only significant bonds threshold distance for PCN costruction : ") or 8.0)      

for protein in proteins_list:
        
    p_name = protein[:4]
    protein_path = proteins_path+p_name+".pdb"
    chains = pcn_final.getAllChainsFromPDB(protein_path)
    chain = None
    if (chains.shape[0]>2):
        chain = str(input("Select the chain to extract from {}: ".format(chains)))[0]
        if chain.upper() not in chains:
            raise Exception("Chain {} not in {}".format(chain, chains))
    index = np.argwhere(chains == chain)
    chains_to_delete = np.delete(chains, index)
    atoms = pcn_final.readPDBFile(protein_path, chain)
    residues = pcn_final.getResidueDistance(atoms, chain)
    dict_residue_name = pcn_final.associateResidueName(residues)
    residue_names = np.array(list (dict_residue_name.items()))
    
    print("computing adjacency matrix for protein {}... (This may take time)".format(p_name))
    A = pcn_final.adjacent_matrix(output_path, residues, p_name, min_, max_)
    
    G = nx.from_numpy_matrix(A)  
    
    for algorithm_choice in algorithms_choice:
               
        method_to_call = getattr(pcn_final, algorithm_choice)
                       
        if (algorithm_choice == 'asyn_fluidc'):                        
            for k in ks:
                labels = method_to_call(G, k)
                n_coms = int( max(labels) + 1)
                pcn_final.save_labels(output_path, labels, residue_names, p_name,  method=algorithm_choice)
                pcn_pymol_scripts.pymol_plot_chain(protein_path, output_path, "Communities", algorithm_choice, n_coms, chains_to_delete)
                       
        else:
            labels = method_to_call(G)
            n_coms = int( max(labels) + 1)
            pcn_final.save_labels(output_path, labels, residue_names, p_name,  method=algorithm_choice)
            pcn_pymol_scripts.pymol_plot_chain(protein_path, output_path, "Communities", algorithm_choice, n_coms, chains_to_delete)