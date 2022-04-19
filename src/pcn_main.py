"""
PCN-Miner: A Software for the analysis of Protein Contact Networks.

Choose the proteins you want to study and insert the PDB codes: for each protein, PCN-Miner will calculate the PCN, then the user can do:
- Node centrality analysis with betweenness, degree, closeness and eigenvector centrality algorithms.
- Structural Analysis of the PCN with Community Detection, Spectral Custering, Embedding+Clustering algorithms.
    1. Supported Spectral Clustering methods: Shi Malik, Normalized an Unnormalized approach both with soft (Fuzzy C-Means) and hard (KMeans) clustering algorithms;
    2. Supported Community Detection methods: Louvain, Leiden, Walktrap, Spinglass, Infomap (not available for Windows users), Async FluidC, Greedy Modularity;
    3. Supported Embeddings methods: HOPE, Laplacian Eigenmaps, Node2Vec.
  
The outputs (node centralities/communities/clusters) will be plotted on the protein structure using the software PyMOL.

Authors: Lomoio Ugo
         Pietro Hiram Guzzi
         Ferrarelli Giuseppe
         
License: CC0-1.0
"""

import os
from networkx import from_numpy_matrix
import pcn_miner
import numpy as np
from scipy.linalg import eigh
import pcn_pymol_scripts
from sys import platform
import configparser

print("Protein Contact Network Miner")

print("Software Available under CC-BY Licence ")

print("Free for Academic Usage")

print("Magna Graecia University of Catanzaro")

print(" ")

input('Press Enter to Continue')

#supported algorithms here 
supported_algorithms_clustering = ["unnorm_ssc", "norm_ssc", "unnorm_hsc", "norm_hsc", "hsc_shimalik", "ssc_shimalik", "skl_spectral_clustering"]
supported_algorithms_embeddings = [
                                   "fuzzycmeans_hope", "kmeans_hope", "fuzzycmeans_laplacianeigenmaps", "kmeans_laplacianeigenmaps" ,
                                   "fuzzycmeans_node2vec", "kmeans_node2vec"
                                  ]
supported_algorithms_communities = ["louvain", "leiden", "walktrap", "asyn_fluidc", "greedy_modularity", "infomap", "spinglass"]

supported_centralities_measures = ["closeness", "eigenvector_c", "betweenness", "degree_c"]

#handle differents slash path separator: example "\" for windows and "/" for Linux and Mac.
#need this for dynamic Output path creation 
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
    print("Infomap algorithm not supported on Windows")
    supported_algorithms_communities.remove("infomap")

#create a dictionary of all supported algoritmhs 
#with this tecnique, the user can choose one, or a list, of supported algorithms without typing the entire name but using the key of the dictionary. 
#the user can choose all the supported algorithms by typing 0. 
"""
dict = {
        0 = "all" 
        1 = ALGORITHM1
        2 = ALGORITHM2
        ...
       }
"""
dict_supported_algorithms_clustering = dict()
dict_supported_algorithms_embeddings = dict()
dict_supported_algorithms_communities = dict()

dict_supported_algorithms_clustering[str (0)] = "all"    #choose all supported algoritmhs 
#populate the dictionaries 
for i, algorithm in enumerate(supported_algorithms_clustering):
    dict_supported_algorithms_clustering[str (i+1)] = algorithm

dict_supported_algorithms_embeddings[str (0)] = "all"
for i, algorithm in enumerate(supported_algorithms_embeddings):
    dict_supported_algorithms_embeddings[str (i+1)] = algorithm

dict_supported_algorithms_communities[str (0)] = "all"   
for i, algorithm in enumerate(supported_algorithms_communities):
    dict_supported_algorithms_communities[str (i+1)] = algorithm

#check and read the config file: in the config file are saved 3 paths: PDB files path, Outputs path, Adjacency matrixs path. 
use_config_choice = None #user can choose to use the paths saved in the config file or to choose 3 new paths
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

while (end==False): #for multiple analysis

    print('As First Step you must choose Format File Input: PDB structures or Preprocessed PCN') #start with PDB files or with a preprocessed PCN?
    initial_choice = str (input("Digit 'pdb' to use .pdb files or  'adj' to load existing PCN: ")).casefold() 
    
    #check if the given output path exists
    if(not os.path.isdir(output_path)):
        os.makedirs(output_path)
    #if not exists, create it 
    if (not output_path.endswith(add_slash_to_path)):
        output_path = output_path+add_slash_to_path
    print('')
    
    #check if the given pdb file path exists 
    is_dir_prot = os.path.isdir(proteins_path)
    
    #if the given pdb file path exists
    if is_dir_prot:
        
        #create a list of all the .pdb files inside the pdb file path   
        if (not proteins_path.endswith(add_slash_to_path)):
            protein_list_dir = [file[:(len(file)-4)] for file in os.listdir(proteins_path) if file.endswith('.pdb')]
            print("List of proteins in {}: {}".format(proteins_path, protein_list_dir))
            
            #choose the pdb files
            print('Please Insert Protein PDB Identifiers, separated by comma, without .pdb, e.g. 7nxc for 7nxc.pdb ')
            protein_choice = str(input("Enter proteins list items (splitted by ',') or 'all': "))
            
            #if the user choose the option "all", he wants to study all the proteins in the pdb file path
            if protein_choice.casefold() == 'all':
                proteins_list = protein_list_dir            
            else:
                proteins_list = [protein.casefold() for protein in protein_choice.replace(" ","").split(",")] #strip space
            
            proteins_path = proteins_path+add_slash_to_path
            pdb_list = proteins_list.copy()
            for i, protein in enumerate(pdb_list):
                if (not protein.endswith('.pdb')):
                    pdb_list[i] = protein+'.pdb'
            #check if the user given pdb codes are in the pdb file path, if they are not here the software will download them
            pcn_miner.checkIfFilesExists(pdb_list, "pdb", proteins_path) 
                    
    else:
         raise Exception("'{}' is not a directory.".format(proteins_path))
    
    #min_ threshold for extracting only non convalent interactions within a.a. of the protein 
    #max_ threshold for extracting only significant interactions within a.a. of the protein 
    min_ = float(input("Entering non covalent bonds threshold distance for PCN costruction: ") or 4.0)    
    max_ = float(input("Entering only significant bonds threshold distance for PCN costruction : ") or 8.0)
    
    #if the user want to analyze an existing PCN
    if(initial_choice == 'adj'):
        
        if(use_config_choice==0):
            print('Please insert the path of the directory containing Adjacency Matrixs')
            adj_filespath = str( input("Insert Adjacency matrix filepath: "))
        is_dir_adj = os.path.isdir(adj_filespath)   #check if the given adj files path exists
        
        #if not exists, create it
        if (not adj_filespath.endswith(add_slash_to_path)):
            adj_filespath = adj_filespath+add_slash_to_path
        
        #if the adj file path and the protein file path exists
        if (is_dir_adj and is_dir_prot):
            adj_list = [protein.casefold()+"_adj_mat_{}_{}.txt".format(min_, max_) for protein in proteins_list]
            pcn_miner.checkIfFilesExists(adj_list, "adj", proteins_path, adj_filespath) #check if really exist the PCN for the given pdb file 
        else:
            raise Exception("'{}' is not a directory.".format(adj_filespath))
            
    elif (initial_choice == 'pdb'):
        pass
    else: 
        raise Exception("'initial_choice' input must be 'pdb' or 'adj' but '{}' given.".format(initial_choice))

    centrality_initial_choice =  int(input("Press 0 if you want to compute a node centrality measure: ")) #ask if the user want to do a node centrality analysis
    
    #structural analysis of the PCN 
    print('Please select the analysis approach: SPECTRAL|EMBEDDINGS|COMMUNITY'+'\n')
    print('Spectral will apply spectral clustering on PCN'+'\n')
    print('Embedding will apply embedding followed by clustering'+'\n')
    print('Community will apply community discovery'+'\n')
    
    #user have to select only one type of analysis
    type_choice = str( input("Choice between 'spectral' (Spectral Clustering), 'embeddings' (Embeddings+Clustering) and 'community' (Community Extraction): ")).casefold()
    
    #create a list of all the selected algorithms 
    algorithms_choice = []
    if (type_choice == 'spectral'):
        dict_supported_algorithms_type_selected = dict_supported_algorithms_clustering
        string_to_print = ""
        for i, algorithm in dict_supported_algorithms_clustering.items():
            string_to_print =  string_to_print+"Type '{}' for '{}' algorithm/s \n".format(i, algorithm)
        algorithms_choice_numeric =  str( input("Choice one or a list (splitted by ',') Spectral Clustering algoritms from: \n{}Choice: ".format(string_to_print)))
        if(algorithms_choice_numeric  == "0"): #if the user type 0, he wants to use all the supported algorithms 
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
    
    #if the user type a number, or a list of integer numbers, not equal to 0, split the string to create a list
    if (algorithms_choice_numeric != "0"):
        if (isinstance(algorithms_choice_numeric, str)):
            if(algorithms_choice_numeric.split(',')):
                algorithms_choice_numeric =  [str(algorithm_choice_numeric) for algorithm_choice_numeric in algorithms_choice_numeric.replace(" ","").split(",") if(algorithm_choice_numeric != "0")]
        
        for algorithm_choice_numeric in algorithms_choice_numeric :
            algorithms_choice.append(dict_supported_algorithms_type_selected[algorithm_choice_numeric])            
        
        print('Selected algorithms: '+ str(algorithms_choice))
    
    #if the user want to use a spectral or embedding clustering approach
    if ((type_choice == 'spectral') or (type_choice == 'embeddings')):
        
        #ask if the user want to use the seme number of clusters for all the proteins 
        k_initial_choice = int(input("Enter 0 if you want to use the same number of clusters k for Spectral Clustering to all the proteins: ") or 1)
        if (k_initial_choice == 0):
            
            #insert the parameters 
            if(type_choice == 'spectral'):
                #number of clusters to extract 
                k_choice = str(input("Entering k for spectral clustering: Enter an int, a list of ints (split with ',') or type 'best_k': ") or 'best_k')  
            else: #embeddings
                #number of clusters 
                k_choice = str(input("Entering k for embedding + clustering: Enter an int, a list of ints (split with ','): "))                                             
                #dimension of embedding
                d = int (input("Enter d parameter for d-dimensional embedding: ") or 2)
                #decay factor, for HOPE embedding
                beta = None
                #random walk lenght, for node2vec
                walk_len = None
                #number of random walks for each node, for node2vec
                num_walks = None
                #for each algorithm in the list of selected algorithms 
                for algorithm_choice in algorithms_choice:
                    #if the algorithm selected is HOPE, insert the beta parameter
                    if ("hope" in algorithm_choice):
                        beta = float(input("Enter beta parameter for d-dimensional HOPE embedding: ") or 0.01)
                        break
                    #if the algorithm selected is node2vec, insert walk_len and num_walks parameters 
                    if("node2vec" in algorithm_choice):
                        walk_len = int(input("Enter the lenght of each random walk: ") or 100)
                        num_walks = int(input("Enter the number of walks per node: ") or 100)
                        break
            
            #if the user want to compute the best number of clusters for spectral clustering, use the max eigengap method
            if ((k_choice == 'best_k') and (type_choice == 'spectral')):
                #number of best k to use
                n_of_best_ks = int(input("Enter the number of best_ks to try: ") or 1)              
            #if the user insert a number of a list of numbers, split the string to create the list of ks       
            elif(k_choice.split(',')):
                ks =  [int(item) for item in k_choice.replace(" ","").split(",")]
                                  
            else:
                raise Exception("'k_choice' input must be an int, a list of ints or 'best_k' but '{}' given.".format(k_choice))   
    #if the user want to do a community extraction analysis        
    elif (type_choice == 'community'):
        
        #if the community detection algorithm selected is asyn_fluidc, he requires a fixed number of communities 
        if ('asyn_fluidc' in algorithms_choice):
            #ask if the user want to use the same number of communities for all the proteins     
            k_initial_choice = int(input("Enter 0 if you want to use the same number of community k for Asyn FluidC to all the proteins: ") or 1)
                                         
            if (k_initial_choice == 0):
                k_choice = str(input("Entering k for Asyn FluidC: Enter an int, a list of ints (split with ','): "))
                #split the string to create a list
                if(k_choice.split(',')):
                    ks =  [int(item) for item in k_choice.replace(" ","").split(",")]
    
    #ask if the user want to compute the partecipation coefficient plot
    if (len(algorithms_choice)>0):
        plot_p = int(input("Press 0 if you want to compute the participation Coef plot: "))
    
    #if the user want to do a centrality analysis
    if (centrality_initial_choice  == 0):
            
            #select the node centrality algorithm/s 
            centralities_choice = str(input("Select one, a list (splitted by ','), or 'all' centrality measures from {}: ".format(str(supported_centralities_measures)))).casefold()     
            #if 'all', select all supported centrality measures
            if centralities_choice == "all":
                centralities_choice = supported_centralities_measures
            #else, try to split the string to create a list of centrality algorithms
            elif(centralities_choice.split(',')):
                centralities_choice =  [str(centrality_choice) for centrality_choice in centralities_choice.replace(" ","").split(",")]
            else:
                raise Exception("{} not supported".format(centralities_choice))
    
    #for each protein in the selected proteins list
    for protein in proteins_list:
        
        p_name = protein
        protein_path = proteins_path+p_name+".pdb"
        atoms = pcn_miner.readPDBFile(protein_path) #read 
        residues = pcn_miner.getResidueCoordinates(atoms) 
        dict_residue_name = pcn_miner.associateResidueName(residues)
        residue_names = np.array(list (dict_residue_name.items()))
        
        #if the user starts with a pdb file, we have to compute the PCN 
        if(initial_choice == 'pdb'):
                     
            print("computing adjacency matrix for protein {}... (This may take time)".format(p_name))
            A = pcn_miner.adjacent_matrix(output_path, residues, p_name, min_, max_)
      
        #if the user starts with a precomputed PCN, read the adj file that represent the PCN
        else:     #'adj'           
            print("reading adjacency matrix for protein {}...".format(p_name))
            A = pcn_miner.read_adj_mat(adj_filespath, p_name, min_, max_)
        
        #if the user wants to do a centrality measure 
        if (centrality_initial_choice  == 0):
            
            print(("protein {}: CENTRALITY MEASURES COMPUTING NOW").format(p_name))
            
            #create the PCN from the adj file 
            G = from_numpy_matrix(A)  
            #extract the residue names: example ALA1
            residue_names_1 = np.array(residue_names[:, 1], dtype = str)        
            
            #for each centrality measure in the selected centrality measures
            for centrality_choice in centralities_choice:

                if(centrality_choice in supported_centralities_measures):
                    #compute the nodes centrality for the graph F
                    print("Computing {} centrality measure on {} PCN".format(centrality_choice, p_name))
                    method_to_call = getattr(pcn_miner, centrality_choice) 
                    centrality_measures = method_to_call(G, residue_names_1)#call the supported method from the pcn_miner file
                    pcn_miner.save_centralities(output_path, centrality_measures, p_name, centrality_choice) #save a txt file 
                    pcn_pymol_scripts.pymol_plot_centralities(output_path, centrality_measures, protein_path, centrality_choice) #plot and save centralities with pymol
                            
                else:
                    print("Centrality method {} not supported".format(centrality_choice))
        
        #if the number of selected algorithm for structural analysis (spectral clusterin, community detection, etc) is greater than 0
        if (len(algorithms_choice)>0):
            
            #compute the PCN from the Adj matrix 
            G = from_numpy_matrix(A)
            #for each algorithm in the selected structural algorithms list
            for algorithm_choice in algorithms_choice:  
                  
                print(("protein {} with algorithm {}: COMPUTING NOW").format(p_name, algorithm_choice))
                
                #if "clustering" (spectral or embeddings) 
                if ((type_choice == 'spectral') or (type_choice == 'embeddings')):  
                    #if the user didn't want to use the same value of k for all the proteins
                    if (k_initial_choice != 0):
                        #ask the number of clusters to use 
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
                    
                    #if the user wants to use the best number of clusters for a spectral clustering, use the max eigengap method 
                    if((k_choice == 'best_k') and (type_choice == 'spectral')):
                        if('shimalik' in algorithm_choice): #if the algorithm follows the Shi Malik approach
                            L = pcn_miner.compute_laplacian_matrix(A)       #unnormalized laplacian matrix
                            D = pcn_miner.degree_matrix(A)                  #degree matrix            
                            eigenvalues = eigh(L, D, eigvals_only=True)     #Shi Malik approach: generalized eigenvalue problem
                                    
                        elif('norm' in algorithm_choice):  #if the algorithm selected follows a normalized approach 
                            L = pcn_miner.compute_normalized_laplacian(A)   #normalized laplacian matrix
                            eigenvalues, eigenvectors  = np.linalg.eig(L)   #compute eigenvectors and eigenvalues of norm L
                                   
                        else: #if the algorithm selected follows an unormalized approach 
                            L = pcn_miner.compute_laplacian_matrix(A)       #unnormalized laplacian matrix
                            eigenvalues, eigenvectors = np.linalg.eig(L)    #compute eigenvectors and eigenvalues of unnorm L
                        
                        #call max eigengap method
                        ks = pcn_miner.computeBestK(eigenvalues, n_k=n_of_best_ks) 
                    
                    print("Selected ks: {}".format(str(ks)))
                    
                    #for each k in the list of number of clusters to try
                    for k in ks:   
                        
                        #call the selected method, with the selected parameter, from the pcn_miner file
                        method_to_call = getattr(pcn_miner, algorithm_choice)
                        if type_choice == 'embeddings':
                            print("{} with {} with k = {}, d = {}, beta = {}, walk_len = {}, num_walks = {}".format(p_name, algorithm_choice, k, d, beta, walk_len, num_walks))
                            labels = method_to_call(A, n_clusters=k, d=d, beta=beta, walk_len=walk_len, num_walks=num_walks)
                        elif type_choice == 'spectral':
                            print("{} with {} with k = {}".format(p_name, algorithm_choice, k))
                            labels = method_to_call(A, n_clusters=k)
                            d=None
                            beta=None
                            walk_len=None
                            num_walks=None
                        #save communities/clusters as a txt file
                        pcn_miner.save_labels(output_path, labels, residue_names, p_name, algorithm_choice, d, beta, walk_len, num_walks)
                       
                        #if the algorithm selected follows a Soft Spectral Clustering approach
                        if "ssc" in algorithm_choice:
                            #not always the real number of clusters extracted with Fuzzy C-Means is equal with the selected number of cluster k 
                            print("Given k = {} but soft clustering algoritmh found k {} clusters".format(k, int(max(labels)+1)))
                            k = int(max(labels)) + 1
                        
                        #pymol plots  
                        if(type_choice == 'embeddings'):
                            pcn_pymol_scripts.pymol_plot_embeddings(protein_path, output_path, "ClustersEmbeddings", algorithm_choice, k, d, beta, walk_len, num_walks)
                                     
                        else:#clustering
                            pcn_pymol_scripts.pymol_plot(protein_path, output_path, "Clusters", algorithm_choice, k)
                        
                        #if the user want to compute the partecipation coefficients
                        if (plot_p == 0): 
                            G = from_numpy_matrix(A) #maybe delete that
                            residue_names_1 = np.array(residue_names[:, 1], dtype = str)
                            p = pcn_miner.participation_coefs(G, labels, residue_names_1)
                            pcn_miner.save_part_coef(output_path, p, p_name, algorithm_choice, k)
                            output_path_p = "{}{}{}{}".format(output_path, add_slash_to_path, algorithm_choice, add_slash_to_path)
                            pcn_pymol_scripts.pymol_plot_part_coefs(p, protein_path, output_path_p, algorithm_choice, k)
                                    
                else:#type_choice = 'community'
                    
                    #call the method 
                    method_to_call = getattr(pcn_miner, algorithm_choice)
                    
                    #if the algorithm is asyn fluidic
                    if (algorithm_choice == 'asyn_fluidc'):
                        #if the user didn't want to use the same number of communities
                        if (k_initial_choice != 0):
                            k_choice = str(input("Entering k for Asyn FluidC: Enter an int, a list of ints (split with ','): "))
                            if(k_choice.split(',')):
                                ks =  [int(item) for item in k_choice.replace(" ","").split(",")]
                        #for each number of communities in the list of numbers of communities to try
                        for k in ks:
                            labels = method_to_call(G, k) #call the method
                            pcn_miner.save_labels(output_path, labels, residue_names, p_name,  method=algorithm_choice) #save the communities as txt file
                            pcn_pymol_scripts.pymol_plot(protein_path, output_path, "Communities", algorithm_choice, k) #plot and save the communities with pymol
                            
                            #if the user wants to compute the partecipation coefficients
                            if (plot_p == 0): 
                                residue_names_1 = np.array(residue_names[:, 1], dtype = str)
                                p = pcn_miner.participation_coefs(G, labels, residue_names_1)
                                pcn_miner.save_part_coef(output_path, p, p_name, algorithm_choice, k) #save the part coefs as txt file
                                output_path_p = "{}{}{}{}".format(output_path, add_slash_to_path, algorithm_choice, add_slash_to_path)
                                pcn_pymol_scripts.pymol_plot_part_coefs(p, protein_path, output_path_p, algorithm_choice, k) #plot and save part coefs with pymol
                    
                    else:#if the community detection algorithm is not Asyn Fluidc, no need to specify the number of communities
                        labels = method_to_call(G) #call the method 
                        n_coms = int( max(labels) + 1)
                        pcn_miner.save_labels(output_path, labels, residue_names, p_name,  method=algorithm_choice) #save communities as txt 
                        pcn_pymol_scripts.pymol_plot(protein_path, output_path, "Communities", algorithm_choice, n_coms) #plot and save communities with pymol
                        
                        #if the user wants to compute the partecipation coefficients
                        if (plot_p == 0): 
                            residue_names_1 = np.array(residue_names[:, 1], dtype = str)
                            p = pcn_miner.participation_coefs(G, labels, residue_names_1)
                            pcn_miner.save_part_coef(output_path, p, p_name, algorithm_choice, n_coms)
                            output_path_p = "{}{}{}{}".format(output_path, add_slash_to_path, algorithm_choice, add_slash_to_path)
                            pcn_pymol_scripts.pymol_plot_part_coefs(p, protein_path, output_path_p, algorithm_choice, n_coms)
                        
                            
    print('Computation Completed.')
    #ask if the user wants to do another analysis (another iteration)
    choice=int(input('Please select 0 to make another analysis, otherwise select another numeric key to end the program: '))
    if (choice!=0):
        end=True
    else:
        proteins_path = proteins_path[:-1]
