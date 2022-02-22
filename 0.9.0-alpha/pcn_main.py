import os
import networkx as nx
import pcn_final
import numpy as np
from scipy.linalg import eigh
import pcn_pymol_scripts
from sys import platform
import configparser


print("Protein Contact Network Miner 0.0.9b ")

print("Software Available under CC-BY Licence ")

print("Free for Academic Usage")

print(" ")

input('Press Enter to Continue')

supported_algorithms_clustering = ["unnorm_ssc", "norm_ssc", "unnorm_hsc", "norm_hsc", "hsc_shimalik", "ssc_shimalik", "skl_spectral_clustering"]
supported_algorithms_embeddings = [
                                   "unnorm_ssc_hope", "norm_ssc_hope", "unnorm_hsc_hope", "norm_hsc_hope", "hsc_shimalik_hope", "ssc_shimalik_hope",
                                   
                                   "unnorm_ssc_laplacianeigenmaps", "norm_ssc_laplacianeigenmaps", "unnorm_hsc_laplacianeigenmaps", "norm_hsc_laplacianeigenmaps",
                                   "hsc_shimalik_laplacianeigenmaps", "ssc_shimalik_laplacianeigenmaps", 
                                   
                                   "unnorm_ssc_node2vec", "norm_ssc_node2vec", "unnorm_hsc_node2vec", "norm_hsc_node2vec", "hsc_shimalik_node2vec", 
                                   "ssc_shimalik_node2vec",
                                   ]
supported_algorithms_communities = ["louvain", "leiden", "walktrap", "asyn_fluidc", "greedy_modularity", "infomap", "spinglass"]

supported_centralities_measures = ["closeness", "eigenvector", "betweenness", "degree_c"]

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


dict_supported_algorithms_clustering = dict()
dict_supported_algorithms_embeddings = dict()
dict_supported_algorithms_communities = dict()

dict_supported_algorithms_clustering[str (0)] = "all"
for i, algorithm in enumerate(supported_algorithms_clustering):
    dict_supported_algorithms_clustering[str (i+1)] = algorithm

dict_supported_algorithms_embeddings[str (0)] = "all"
for i, algorithm in enumerate(supported_algorithms_embeddings):
    dict_supported_algorithms_embeddings[str (i+1)] = algorithm

dict_supported_algorithms_communities[str (0)] = "all"   
for i, algorithm in enumerate(supported_algorithms_communities):
    dict_supported_algorithms_communities[str (i+1)] = algorithm

use_config_choice = None
if(os.path.isfile(os.getcwd()+add_slash_to_path+"config.ini")):
    print("config file found!")
    config_obj = configparser.ConfigParser()
    config_obj.read(os.getcwd()+add_slash_to_path+"config.ini")
    paths = config_obj["user_paths"]
    output_path = paths["output_path"]
    proteins_path = paths["proteins_path"]
    adj_filespath = paths["adj_filespath"]
    
    print("Paths in the config file: \n output path = {} \n proteins_path = {} \n adj_filespath = {}".format(output_path, proteins_path, adj_filespath))
    use_config_choice=int(input("Please select 0 if you DON'T want to use the paths in this config file, otherwise select another numeric key: "))
    
    if(use_config_choice==0):
            
        print('Input the Directory in which you want to store the outputs')
        print('The software will create three subdirectories')
        output_path = str (input("Insert Root Path of the Outputs: "))
        print('Input the Directory containing Input Files') 
        proteins_path = str (input("Insert Proteins filepath: "))

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
        
    print('config file successfully created')

end=False 

while (end==False):

    print('As First Step you must choose Format File Input: PDB structures or Preprocessed PCN')
    initial_choice = str (input("Digit 'pdb' to use .pdb files or  'adj' to load existing PCN: ")).casefold()
    
    if(not os.path.isdir(output_path)):
        os.makedirs(output_path)
    
    if (not output_path.endswith(add_slash_to_path)):
        output_path = output_path+add_slash_to_path
    print('')
    
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
                    
    else:
         raise Exception("'{}' is not a directory.".format(proteins_path))
    
    min_ = float(input("Entering non covalent bonds threshold distance for PCN costruction: ") or 4.0)    
    max_ = float(input("Entering only significant bonds threshold distance for PCN costruction : ") or 8.0)
    
    if(initial_choice == 'adj'):
        
        if(use_config_choice==0):
            print('Please insert the path of the directory containing Adjacency Matrixs')
            adj_filespath = str( input("Insert Adjacency matrix filepath: "))
        is_dir_adj = os.path.isdir(adj_filespath)
        
        if (not adj_filespath.endswith(add_slash_to_path)):
            adj_filespath = adj_filespath+add_slash_to_path
            
        if (is_dir_adj and is_dir_prot):
            adj_list = [protein.casefold()+"_adj_mat_{}_{}.txt".format(min_, max_) for protein in proteins_list]
            pcn_final.checkIfFilesExists(adj_list, "adj", proteins_path, adj_filespath) #add_slash_to_path is path separator 
        else:
            raise Exception("'{}' is not a directory.".format(adj_filespath))
            
    elif (initial_choice == 'pdb'):
        pass
    else: 
        raise Exception("'initial_choice' input must be 'pdb' or 'adj' but '{}' given.".format(initial_choice))

    centrality_initial_choice =  int(input("Press 0 if you want to compute a node centrality measure: "))
    
    print('Please select the analysis approach: SPECTRAL|EMBEDDINGS|COMMUNITY'+'\n')
    print('Spectral will apply spectral clustering on PCN'+'\n')
    print('Embedding will apply embedding followed by clustering'+'\n')
    print('Community will apply community discovery'+'\n')
    
    type_choice = str( input("Choice between 'spectral' (Spectral Clustering), 'embeddings' (Embeddings+Clustering) and 'community' (Community Extraction): ")).casefold()
    algorithms_choice = []
    if (type_choice == 'spectral'):
        dict_supported_algorithms_type_selected = dict_supported_algorithms_clustering
        string_to_print = ""
        for i, algorithm in dict_supported_algorithms_clustering.items():
            string_to_print =  string_to_print+"Type '{}' for '{}' algorithm/s \n".format(i, algorithm)
        algorithms_choice_numeric =  str( input("Choice one or a list (splitted by ',') Spectral Clustering algoritms from: \n{}Choice: ".format(string_to_print)))
        if(algorithms_choice_numeric  == "0"): #'all'
            algorithms_choice = supported_algorithms_clustering
            print('Selected algorithms: '+ str(algorithms_choice))
        
                
    elif (type_choice == 'embeddings'):
        dict_supported_algorithms_type_selected = dict_supported_algorithms_embeddings
        string_to_print = ""
        for i, algorithm in dict_supported_algorithms_embeddings.items():
            string_to_print =  string_to_print+"Type '{}' for '{}' algorithm/s \n".format(i, algorithm)
        algorithms_choice_numeric =  str( input("Choice one or a list (splitted by ',') Clustering Embeddings algoritms from: \n{}Choice: ".format(string_to_print)))
        if(algorithms_choice_numeric == "0"):
            algorithms_choice =  supported_algorithms_embeddings
            print('Selected algorithms: '+ str(algorithms_choice))

    elif (type_choice == 'community'):
        dict_supported_algorithms_type_selected = dict_supported_algorithms_communities
        string_to_print = ""
        for i, algorithm in dict_supported_algorithms_communities.items():
            string_to_print =  string_to_print+"Type '{}' for '{}' algorithm/s \n".format(i, algorithm)
        algorithms_choice_numeric = str( input("Choice one or a list (splitted by ',') Community Extraction algoritms from: \n{}Choice: ".format(string_to_print)))
        if(algorithms_choice_numeric == "0"):
            algorithms_choice = supported_algorithms_communities
            print('Selected algorithms: '+ str(algorithms_choice))
    else:
        raise Exception("'type_choice' input must be 'spectral', 'embeddings' or 'community' but '{}' given.".format(type_choice))           
    
    if (algorithms_choice_numeric != "0"):
        if (isinstance(algorithms_choice_numeric, str)):
            if(algorithms_choice_numeric.split(',')):
                algorithms_choice_numeric =  [str(algorithm_choice_numeric) for algorithm_choice_numeric in algorithms_choice_numeric.replace(" ","").split(",") if(algorithm_choice_numeric != "0")]
        
        for algorithm_choice_numeric in algorithms_choice_numeric :
            algorithms_choice.append(dict_supported_algorithms_type_selected[algorithm_choice_numeric])            
        
        print('Selected algorithms: '+ str(algorithms_choice))
                   
    if ((type_choice == 'spectral') or (type_choice == 'embeddings')):
           
        k_initial_choice = int(input("Enter 0 if you want to use the same number of clusters k for Spectral Clustering to all the proteins: ") or 1)
        if (k_initial_choice == 0):
              
            if(type_choice == 'spectral'):
                k_choice = str(input("Entering k for spectral clustering: Enter an int, a list of ints (split with ',') or type 'best_k': ") or 'best_k')  
            else: #embeddings
                k_choice = str(input("Entering k for embedding + clustering: Enter an int, a list of ints (split with ','): "))                                             
                d = int (input("Enter d parameter for d-dimensional embedding: ") or 2)
                beta = None
                walk_len = None
                num_walks = None
                for algorithm_choice in algorithms_choice:
                    if ("hope" in algorithm_choice):
                        beta = float(input("Enter beta parameter for d-dimensional HOPE embedding: ") or 0.01)
                        break
                    if("node2vec" in algorithm_choice):
                        walk_len = int(input("Enter the lenght of each random walk: ") or 100)
                        num_walks = int(input("Enter the number of walks per node: ") or 100)
                        break
              
            if ((k_choice == 'best_k') and (type_choice == 'spectral')):
                n_of_best_ks = int(input("Enter the number of best_ks to try: ") or 1)              
                   
            elif(k_choice.split(',')):
                ks =  [int(item) for item in k_choice.replace(" ","").split(",")]
                                  
            else:
                raise Exception("'k_choice' input must be an int, a list of ints or 'best_k' but '{}' given.".format(k_choice))   
            
    elif (type_choice == 'community'):
            
        if ('asyn_fluidc' in algorithms_choice):
                
            k_initial_choice = int(input("Enter 0 if you want to use the same number of community k for Asyn FluidC to all the proteins: ") or 1)
                                         
            if (k_initial_choice == 0):
                k_choice = str(input("Entering k for Asyn FluidC: Enter an int, a list of ints (split with ','): "))
            
                if(k_choice.split(',')):
                    ks =  [int(item) for item in k_choice.replace(" ","").split(",")]
   
    if (len(algorithms_choice)>0):
        plot_p = int(input("Press 0 if you want to compute the participation Coef plot: "))
    
    for protein in proteins_list:
        
        p_name = protein[:4]
        protein_path = proteins_path+p_name+".pdb"
        atoms = pcn_final.readPDBFile(protein_path)
        residues = pcn_final.getResidueDistance(atoms)
        dict_residue_name = pcn_final.associateResidueName(residues)
        residue_names = np.array(list (dict_residue_name.items()))
        
        if(initial_choice == 'pdb'):
                     
            print("computing adjacency matrix for protein {}... (This may take time)".format(p_name))
            A = pcn_final.adjacent_matrix(output_path, residues, p_name, min_, max_)
      
        else:     #'adj'           
            print("reading adjacency matrix for protein {}...".format(p_name))
            A = pcn_final.read_adj_mat(adj_filespath, p_name, min_, max_)
            
        if (centrality_initial_choice  == 0):
            
            print(("protein {}: CENTRALITY MEASURES COMPUTING NOW").format(p_name))
                              
            G = nx.from_numpy_matrix(A)  
            centralities_choice = str(input("Select one, a list (splitted by ','), or 'all' centrality measures from {}: ".format(str(supported_centralities_measures)))).casefold()
            residue_names_1 = np.array(residue_names[:, 1], dtype = str)        
            if centralities_choice == "all":
                centralities_choice = supported_centralities_measures
            elif(centralities_choice.split(',')):
                centralities_choice =  [str(centrality_choice) for centrality_choice in centralities_choice.replace(" ","").split(",")]
            else:
                raise Exception("{} not supported".format(centralities_choice))
                                                 
            for centrality_choice in centralities_choice:
            
                print("Computing {} centrality measure on {} PCN".format(centrality_choice, p_name))
                if(centrality_choice in supported_centralities_measures):
                            
                    method_to_call = getattr(pcn_final, centrality_choice)
                                
                    centrality_measures = method_to_call(G, p_name, residue_names_1, 10)
                    print(centrality_measures)
                    pcn_pymol_scripts.pymol_plot_centralities(centrality_measures, protein_path, output_path, centrality_choice)
                            
                else:
                    print("{} not supported".format(centrality_choice))
        
        if (len(algorithms_choice)>0):

            G = nx.from_numpy_matrix(A)
            
            for algorithm_choice in algorithms_choice:  
                  
                print(("protein {} with algorithm {}: COMPUTING NOW").format(p_name, algorithm_choice))
                
                if ((type_choice == 'spectral') or (type_choice == 'embeddings')):  
                
                    if (k_initial_choice != 0):
                        
                        if (type_choice == 'spectral'):
                            k_choice = str(input("Entering k for spectral clustering {} algorithm: Enter an int, a list of ints (split with ',') or type 'best_k': ".format(algorithm_choice)))                               
                            if (k_choice == 'best_k'):     
                                n_of_best_ks = int(input("Enter the number of best_ks to try: "))   
                        else:
                            k_choice = str(input("Entering k for embedding+clustering: Enter an int, a list of ints (split with ','): "))                                            
                        
                        if (k_choice == 'best_k'): 
                            pass
                        elif(k_choice.split(',')):
                            ks =  [int(item) for item in k_choice.replace(" ","").split(",")]
                               
                        else:
                            raise Exception("'k_choice' input must be an int, a list of ints or 'best_k' but '{}' given.".format(k_choice))   
                          
                    if((k_choice == 'best_k') and (type_choice == 'spectral')):
                        if('shimalik' in algorithm_choice):
                            L = pcn_final.compute_laplacian_matrix(A)
                            D = pcn_final.degree_matrix(A) 
                            eigenvalues = eigh(L, D, eigvals_only=True)   
                                    
                        elif('norm' in algorithm_choice):
                            L = pcn_final.compute_normalized_laplacian(A)
                            eigenvalues, eigenvectors  = np.linalg.eig(L)    
                                   
                        else: #'unnorm' in algorithm_choice
                            L = pcn_final.compute_laplacian_matrix(A)
                            eigenvalues, eigenvectors = np.linalg.eig(L)
                                    
                        ks = pcn_final.computeBestK(eigenvalues, n_k=n_of_best_ks) 
                    
                    print("Selected ks: {}".format(str(ks)))
                    
                    for k in ks:   
                        
                        print("{} with {} with k = {}".format(p_name, algorithm_choice, k))
                        method_to_call = getattr(pcn_final, algorithm_choice)
                        if type_choice == 'embeddings':
                            labels = method_to_call(A, n_clusters=k, d=d, beta=beta, walk_len=walk_len, num_walks=num_walks)
                        elif type_choice == 'spectral':
                            labels = method_to_call(A, n_clusters=k)
                            d=None
                            beta=None
                            walk_len=None
                            num_walks=None
                        pcn_final.save_labels(output_path, labels, residue_names, p_name, algorithm_choice, d, beta, walk_len, num_walks)
                        
                        if "ssc" in algorithm_choice:
                            print("Given k = {} but soft clustering algoritmh found k {} clusters".format(k, int(max(labels)+1)))
                            k = int(max(labels)) + 1
                        
                        #pymol 
                        if(type_choice == 'embeddings'):
                            pcn_pymol_scripts.pymol_plot_embeddings(protein_path, output_path, "ClustersEmbeddings", algorithm_choice, k, d, beta, walk_len, num_walks)
                                     
                        else:#clustering
                            pcn_pymol_scripts.pymol_plot(protein_path, output_path, "Clusters", algorithm_choice, k)
                            
                        if (plot_p == 0): 
                            G = nx.from_numpy_matrix(A) 
                            residue_names_1 = np.array(residue_names[:, 1], dtype = str)
                            output_path_p = "{}{}{}{}".format(output_path, add_slash_to_path, algorithm_choice, add_slash_to_path)
                            p = pcn_final.participation_coefs(G, labels, residue_names_1)
                            pcn_pymol_scripts.pymol_plot_part_coefs(p, protein_path, output_path_p, algorithm_choice, k)
                                    
                else:#type_choice = 'community'
                     
                    method_to_call = getattr(pcn_final, algorithm_choice)
                    
                    if (algorithm_choice == 'asyn_fluidc'):
                        if (k_initial_choice != 0):
                            k_choice = str(input("Entering k for Asyn FluidC: Enter an int, a list of ints (split with ','): "))
                            if(k_choice.split(',')):
                                ks =  [int(item) for item in k_choice.replace(" ","").split(",")]
                        
                        for k in ks:
                            labels = method_to_call(G, k)
                            n_coms = int( max(labels) + 1)
                            pcn_final.save_labels(output_path, labels, residue_names, p_name,  method=algorithm_choice)
                            pcn_pymol_scripts.pymol_plot(protein_path, output_path, "Communities", algorithm_choice, n_coms)
                            
                            if (plot_p == 0): 
                                residue_names_1 = np.array(residue_names[:, 1], dtype = str)
                                p = pcn_final.participation_coefs(G, labels, residue_names_1)
                                output_path_p = "{}{}{}{}".format(output_path, add_slash_to_path, algorithm_choice, add_slash_to_path)
                                pcn_pymol_scripts.pymol_plot_part_coefs(p, protein_path, output_path_p, algorithm_choice, k) 
                    
                    else:
                        labels = method_to_call(G)
                        n_coms = int( max(labels) + 1)
                        pcn_final.save_labels(output_path, labels, residue_names, p_name,  method=algorithm_choice)
                        pcn_pymol_scripts.pymol_plot(protein_path, output_path, "Communities", algorithm_choice, n_coms)    
                        
                        if (plot_p == 0): 
                            residue_names_1 = np.array(residue_names[:, 1], dtype = str)
                            p = pcn_final.participation_coefs(G, labels, residue_names_1)
                            output_path_p = "{}{}{}{}".format(output_path, add_slash_to_path, algorithm_choice, add_slash_to_path)
                            pcn_pymol_scripts.pymol_plot_part_coefs(p, protein_path, output_path_p, algorithm_choice, n_coms)
                        
                            
    print('Computation Completed.')
    choice=int(input('Please select 0 to make another analsys, otherwise select another numeric key to end the program: '))
    if (choice!=0):
        end=True
    else:
        proteins_path = proteins_path[:-1]
