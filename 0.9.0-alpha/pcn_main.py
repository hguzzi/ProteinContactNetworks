import os
import networkx as nx
import pcn_final
import numpy as np
from scipy.linalg import eigh
import pcn_pymol_scripts
from sys import platform

supported_algorithms_clustering = ["Unnorm_SSC", "Norm_SSC", "Unnorm_HSC", "Norm_HSC", "HSC_ShiMalik", "SSC_ShiMalik"]
supported_algorithms_embeddings = ["Unnorm_SSC_HOPE", "Norm_SSC_HOPE", "Unnorm_HSC_HOPE", "Norm_HSC_HOPE", "HSC_ShiMalik_HOPE", "SSC_ShiMalik_HOPE",
                                    "Unnorm_SSC_LaplacianEigenmaps", "Norm_SSC_LaplacianEigenmaps", "Unnorm_HSC_LaplacianEigenmaps", "Norm_HSC_LaplacianEigenmaps",
                                    "HSC_ShiMalik_LaplacianEigenmaps", "SSC_ShiMalik_LaplacianEigenmaps"]
supported_algorithms_communities = ["Louvain", "Leiden", "Walktrap", "Asyn_fluidc", "Greedy_modularity", "Infomap", "Spinglass"]


supported_algorithms_clustering_lw = ["unnorm_ssc", "norm_ssc", "unnorm_hsc", "norm_hsc", "hsc_shimalik", "ssc_shimalik"]
supported_algorithms_embeddings_lw=list()
for i in supported_algorithms_embeddings :
       supported_algorithms_embeddings_lw.append(i.lower()) 
supported_algorithms_communities_lw = ["louvain", "leiden", "walktrap", "asyn_fluidc", "greedy_modularity", "infomap", "spinglass"]



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
    supported_algorithms_communities.remove("Infomap")


print("Protein Contact Network Analyzer 0.0.9b ")

print("Software Available under CC-BY Licence ")

print("Free for Academic Usage")

print(" ")

print('Press Enter to Continue')

end=False 

while (end==False):
    print('As First Step you must choose Format File Input: PDB structures or Preprocessed PCN')
    initial_choice = str (input("Digit pdb' to use .pdb files or  'adj' to load existing PCN: "))
    
    
    print('Input the Directory in which you want to store the outputs')
    print('The software will create three subdirectories')
    output_path = str (input("Insert Root Path of the Outputs: "))
    
    if(not os.path.isdir(output_path)):
        os.makedirs(output_path)
    
    if (not output_path.endswith(add_slash_to_path)):
        output_path = output_path+add_slash_to_path
    print('')
    
    print('Input the Directory containing Input Files')
    
    proteins_path = str (input("Insert Proteins filepath: "))
    is_dir_prot = os.path.isdir(proteins_path)
    
    if is_dir_prot:
            
        if (not proteins_path.endswith(add_slash_to_path)):
            proteins_path = proteins_path+add_slash_to_path
            print('Please Insert Protein PDB Identifiers, separated by comma, without .pdb, e.g. 7nxc for 7nxc.pdb ')
            proteins_list = [protein for protein in input("Enter proteins list items (splitted by ','): ").split(',')]
            pdb_list = proteins_list.copy()
            for i, protein in enumerate(pdb_list):
                if (not protein.endswith('.pdb')):
                    pdb_list[i] = protein+'.pdb'
            all_prot_files_exists = pcn_final.checkIfFilesExists(pdb_list, "pdb", proteins_path) #add_slash_to_path is path separator 
                    
    else:
         raise Exception("'{}' is not a directory.".format(proteins_path))
      
    if(initial_choice == 'adj'):
        print('Please insert the path of the directory containing Adjacency Matrices')
        adj_filespath = str( input("Insert Adjacency matrix filepath: "))
    
        is_dir_adj = os.path.isdir(adj_filespath)
        
        if (not adj_filespath.endswith(add_slash_to_path)):
            adj_filespath = adj_filespath+add_slash_to_path
            
        if (is_dir_adj and is_dir_prot):
            adj_list = [protein+"_adj_mat.txt" for protein in proteins_list]
            all_adj_files_exists = pcn_final.checkIfFilesExists(adj_list, initial_choice, adj_filespath)
            if (all_adj_files_exists and all_prot_files_exists):
                all_files_exists=True
            else: 
                all_files_exists = False
        else:
            raise Exception("'{}' is not a directory.".format(adj_filespath))
    elif (initial_choice == 'pdb'):
        all_files_exists = all_prot_files_exists
    else: 
        raise Exception("'initial_choice' input must be 'pdb' or 'adj' but '{}' given.".format(initial_choice))
       
    if all_files_exists:
        print('Please select the analysis approach: SPECTRAL|EMBEDDINGS|COMMUNITY'+'\n')
        print('Spectral will apply spectral clustering on PCN'+'\n')
        print('Embedding will apply embedding followed by clustering'+'\n')
        print('Community will apply community discovery'+'\n')
        type_choice = str( input("Choice between 'spectral' (Spectral Clustering), 'embeddings' (Embeddings+Spectral Clustering) and 'community' (Community Extraction): "))
        if (type_choice.casefold() == 'spectral'):
            algorithm_choice = str( input("Choice one Spectral Clustering algoritm from {}: ".format(str(supported_algorithms_clustering))))
        elif (type_choice.casefold() == 'embeddings'):
            algorithm_choice = str( input("Choice one Embedded Spectral Clustering algoritm from {}: ".format(str(supported_algorithms_embeddings))))
            d = int (input("Enter d parameter for d-dimensional embedding: "))
            beta = None
            if ("HOPE" in algorithm_choice):
                beta = float(input("Enter beta parameter for d-dimensional HOPE embedding: "))
        elif (type_choice == 'community'):
            algorithm_choice = str( input("Choice one Community Extraction algoritm from {}: ".format(str(supported_algorithms_communities))))
            print('Algorithm:'+algorithm_choice)
        else:
            raise Exception("'type_choice' input must be 'spectral', 'embeddings' or 'community' but '{}' given.".format(type_choice))
        
        if ((algorithm_choice not in supported_algorithms_clustering) and (algorithm_choice not in supported_algorithms_embeddings) and
            (algorithm_choice not in supported_algorithms_communities)):
            raise Exception("Algorithm {} not supported yet.".format(algorithm_choice))
        else:
            for protein in proteins_list:
    
                p_name = protein[:4]   #p = 6vxx, protein = 6vxx.pdb, adj = 6vxx_adj_mat.txt
                print(("protein {}: COMPUTING NOW").format(p_name))
                
                protein_path = proteins_path+p_name+".pdb"
                atoms = pcn_final.readPDBFile(protein_path)
                residues = pcn_final.getResidueDistance(atoms)
                dict_residue_name = pcn_final.associateResidueName(residues)
                residue_names = np.array(list (dict_residue_name.items()))
                
                if(initial_choice == 'pdb'):
                    
                    min_ = int(input("Entering non covalent bonds threshold distance for PCN costruction: ") or 4)    
                    max_ = int(input("Entering only significant bonds threshold distance for PCN costruction : ") or 8)
                    print("computing adjacency matrix... (This may take time)")
                    A = pcn_final.adjacent_matrix(output_path, residues, p_name, min_, max_)
    
                else:     #'adj'           
                    print("reading adjacency matrix...")
                    A = pcn_final.read_adj_mat(adj_filespath, p_name)
                    
                
                if ((algorithm_choice in supported_algorithms_clustering) or (algorithm_choice in supported_algorithms_embeddings)):
    
                    k_choice = str(input("Entering k for clustering: Enter an int, a list of ints (split with ',') or type 'best_k': "))    
                    if (k_choice == 'best_k'):
                        
                        n_of_best_ks = int(input("Enter the number of best_ks to try: "))
                    
                        if('ShiMalik' in algorithm_choice):
                            L = pcn_final.compute_laplacian_matrix(A)
                            D = pcn_final.degree_matrix(A) 
                            eigenvalues = eigh(L, D, eigvals_only=True)   
                        
                        elif('Norm' in algorithm_choice):
                            L = pcn_final.compute_normalized_laplacian(A)
                            eigenvalues, eigenvectors  = np.linalg.eig(L)    
                        
                        else: #'Unnorm' in algorithm_choice
                            L = pcn_final.compute_laplacian_matrix(A)
                            eigenvalues, eigenvectors = np.linalg.eig(L)
                        
                        ks = pcn_final.computeBestK(eigenvalues, n_k=n_of_best_ks) 
                       
                    elif(k_choice.split(',')):
                        ks =  [int(item) for item in k_choice.split(',')]
                    
                    else:
                        raise Exception("'k_choice' input must be an int, a list of ints or 'best_k' but '{}' given.".format(k_choice))   
    
                    for k in ks:     
                        method_to_call = getattr(pcn_final, algorithm_choice)
                        
                        if(algorithm_choice in supported_algorithms_embeddings):
                            labels = method_to_call(A, k, d, beta)
                        else: #algorithm_choice in supported_algorithms_clustering
                            labels = method_to_call(A, k)
                       ## /d=None
                       ##     beta=None
    
                        pcn_final.save_labels(output_path, labels, residue_names, p_name, algorithm_choice, d, beta)
                        if(type_choice == 'embeddings'):
                            pcn_pymol_scripts.pymol_plot_embeddings(protein_path, output_path, "ClustersEmbeddings", algorithm_choice, k, d, beta)
                        
                        else:
                            pcn_pymol_scripts.pymol_plot(protein_path, output_path, "Clusters", algorithm_choice, k)
                            
                else:# is community extraction
                    G = nx.from_numpy_matrix(A)
                    method_to_call = getattr(pcn_final, algorithm_choice)
                    labels = method_to_call(G)
                    pcn_final.save_labels(output_path, labels, residue_names, p_name,  method=algorithm_choice)
                    n_coms = int( max(labels) + 1)
                    pcn_pymol_scripts.pymol_plot(protein_path, output_path, "Communities", algorithm_choice, n_coms)
    print('Computation Completed')
    choice=int(input('Please select 0 to make another analsys, 1 to end program'))
    if (choice==0):
         end=False