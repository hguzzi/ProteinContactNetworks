#IMPORT LIBRARIES

import sys
import numpy as np
import ast
import networkx as nx
import urllib
from scipy.spatial.distance import euclidean
from sklearn.cluster import SpectralClustering, KMeans
from fcmeans import FCM
import collections
from operator import itemgetter
from sklearn import metrics
import scipy
from scipy.linalg import eigh
import warnings
import os
import base64 
import matplotlib.pyplot as plt# NEW
#DELETED PRODY LIBRARY

from networkx.algorithms.centrality import degree_centrality, eigenvector_centrality, closeness_centrality, betweenness_centrality

from gem.embedding.hope     import HOPE
from gem.embedding.lap      import LaplacianEigenmaps
#from gem.embedding.node2vec import node2vec

from cdlib import algorithms
#from networkx.algorithms.community.centrality import girvan_newman as girvan_newman_
from networkx.algorithms.community.asyn_fluid import asyn_fluidc as asyn_fluidc_
from networkx.algorithms.community import greedy_modularity_communities

#CREATE Protein Contact Network
from sys import platform

if platform == "linux" or platform == "linux2":
    # linux
    add_slash_to_path = '/'
elif platform == "darwin":
    # OS X
    add_slash_to_path = '/'
elif platform == "win32":
    # Windows...
    add_slash_to_path = '\\' 
       
def checkIfFilesExists(files, initial_choice, proteins_path, adj_path = None):
    
    not_existing_pdb_files = []
    all_pdb_files_exists = True
    
    for file in files:
        p_name = file[:4]
        pdb_path = "{}{}.pdb".format(proteins_path, p_name)
        if((not os.path.isfile(pdb_path))):
            all_pdb_files_exists = False
            not_existing_pdb_files.append(file)
            
    if(not all_pdb_files_exists):
        for protein in not_existing_pdb_files:
            p_name = protein[:4]
            print(("protein {} pdb file missing, fetching on PDB database...").format(p_name))
            urllib.request.urlretrieve("http://files.rcsb.org/download/{}.pdb".format(p_name), "{}{}.pdb".format(proteins_path, p_name))
            print(("protein {} fetched successfully").format(p_name))
        all_pdb_file_exists = True
    
    if(initial_choice == "adj"):
        not_existing_adj_files = []
        all_adj_files_exists = True
        for file in files:
            
            file_path = "{}{}".format(adj_path, file)
            if((not os.path.isfile(file_path))):
                all_adj_files_exists = False
                not_existing_adj_files.append(file)
           
        if(not all_adj_files_exists):
            for protein in not_existing_adj_files:
                p_name = protein[:4]   #adj = 6vxx_adj_mat.txt
                print(("protein {} adj matrix missing...").format(p_name))
                protein_path = proteins_path+p_name+".pdb"
                atoms = readPDBFile(protein_path)
                residues = getResidueDistance(atoms)
                dict_residue_name = associateResidueName(residues)
                residue_names = np.array(list (dict_residue_name.items()))
                                                
                min_ = int(input("Entering non covalent bonds threshold distance for PCN costruction: ") or 4)    
                max_ = int(input("Entering only significant bonds threshold distance for PCN costruction : ") or 8)
                print("computing adjacency matrix... (This may take time)")
                output_path = os.path.abspath(os.path.join(adj_path, os.pardir))+add_slash_to_path
                A = adjacent_matrix(output_path, residues, p_name, min_, max_)
            
            all_adj_file_exists = True                   
    
def readPDBFile(pdbFilePath, chain = None):

    atoms = []
    with open(pdbFilePath) as pdbfile:
        for line in pdbfile:
            if line[:4] == 'ATOM':
              # Split the line
                splitted_line = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54], line[56:61], line[62:66]]
                if chain is None:
                    atoms.append(splitted_line)
                else:
                    if(line[21] == chain.upper()):
                        atoms.append(splitted_line)
                # To format again the pdb file with the fields extracted
                #print("%-6s%5s %4s %3s %s %4s    %8s%8s%8s   %3s%3s\n"%tuple(splitted_line))
    return np.array(atoms)    

def getAllChainsFromPDB(pdbFilePath):

    chains = []
    with open(pdbFilePath) as pdbfile:
        for line in pdbfile:
            if line[:4] == 'ATOM':
              # Split the line
                splitted_line = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54], line[56:61], line[62:66]]
                if(line[21] not in chains):
                    chains.append(line[21])
                # To format again the pdb file with the fields extracted
                #print("%-6s%5s %4s %3s %s %4s    %8s%8s%8s   %3s%3s\n"%tuple(splitted_line))
    return np.array(chains)
    
def getResidueDistance(atoms, chain = None):
  
    residues = []
    residues_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 
                   'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 
                   'TRP', 'TYR']
  
    dropped = []

    last_residue_num = 0

    for i, atom in enumerate(atoms):
        
        residue_name = atom[3]
        residue_chain = atom[4]
        residue_num = int (atom[5])
        
        if residue_name in residues_list:   

            if (residue_num != last_residue_num):
                
                if (atom[2].replace(" ","")=="CA"):#'C'or 'CA'
                    if chain is None:
                        cord_C = atom[6:9] 
                        residues.append([residue_name + str (residue_num) + " " + residue_chain, cord_C]) 
                        last_residue_num = residue_num   
                    else: 
                        if residue_chain == chain:
                            cord_C = atom[6:9] 
                            residues.append([residue_name + str (residue_num) + " " + residue_chain, cord_C]) 
                            last_residue_num = residue_num  
        else:
          dropped.append(residue_name)

    return np.array(residues)
     
def associateResidueName(residues):
  
    dict_residue_name = dict()
    for i in range(residues.shape[0]):
        dict_residue_name[str (i)] = residues[i, 0]
    return dict_residue_name
 
def getResiduesSequence(pbdFilePath, chain = None):
  
    seq_res = []
    residues_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
    
    with open(pbdFilePath) as pdbfile:
        for line in pdbfile:
            if line[:6]=='SEQRES':
                if chain is None:
                    splitted_line = [line[19:22], line[23:26], line[27:30], line[31:34], line[35:38],
                                    line[39:42], line[43:46], line[47:50], line[51:54], line[55:58],
                                    line[59:62], line[63:66], line[67:70]]
                    for residue in splitted_line:
                        if residue in residues_list:
                            seq_res.append(residue)
                else: 
                    if chain == line[12]:
                        if chain is None:
                            splitted_line = [line[19:22], line[23:26], line[27:30], line[31:34], line[35:38],
                                            line[39:42], line[43:46], line[47:50], line[51:54], line[55:58],
                                            line[59:62], line[63:66], line[67:70]]
                            for residue in splitted_line:
                                if residue in residues_list:
                                    seq_res.append(residue)
                         
    return np.array(seq_res)

def read_adj_mat(adj_filepath, p, min_, max_):

    if (os.path.isfile("{}{}_adj_mat_{}_{}.txt".format(adj_filepath, p, min_, max_))):
        adj = np.loadtxt("{}{}_adj_mat_{}_{}.txt".format(adj_filepath, p, min_, max_))
        return adj
    else: 
        raise Exception("Adj matrix for protein {} doesn't exists.".format(p))

#compute adj matrix 
def adjacent_matrix(output_path, residues, p, min_=4, max_=8):
    
    n = residues.shape[0]
    adj = np.zeros((n,n))
    d = np.zeros((n,n), dtype=float)
    edge_list = []

    for i in range(n):
        for j in range(n):     
            if (i!=j): 
                p1 = np.array(residues[i, 1], dtype=float) #CORD_C 
                p2 = np.array(residues[j, 1], dtype=float) #CORD_C
                d[i][j] = euclidean(p1, p2)
                if ((d[i][j]>min_) and (d[i][j]<max_)):
                    adj[i][j] = 1
                    adj[j][i] = 1 #for symmetricity
                    if (([j, i] not in edge_list) and ([i, j] not in edge_list)):
                        edge_list.append([i, j])
    
    
    if not os.path.exists("{}Distances".format(output_path)):
        os.makedirs("{}Distances".format(output_path))
    np.savetxt("{}Distances{}{}_distances.txt".format(output_path, add_slash_to_path, p), d, fmt='%.2f')
    print("saved distances matrix")
    
        
    if not os.path.exists("{}Edgelists".format(output_path)):
        os.makedirs("{}Edgelists".format(output_path))
        
    np.savetxt("{}Edgelists{}{}_edgelist_{}_{}.csv".format(output_path, add_slash_to_path, p, min_, max_), np.array(edge_list), fmt='%.2f')
    print("saved edge list")
        
    if not os.path.exists("{}Adj".format(output_path)):
        os.makedirs("{}Adj".format(output_path))
    np.savetxt("{}Adj{}{}_adj_mat_{}_{}.txt".format(output_path, add_slash_to_path, p, min_, max_), adj, fmt='%.2f')
    print("saved adj matrix")
    
    return adj

#SAVE LABELS
def save_labels(output_path, labels, residue_names, p_name, method=None, d=None, beta=None):

    supported_methods_clustering = ["unnorm_ssc", "norm_ssc", "unnorm_hsc", "norm_hsc", "hsc_shimalik", "ssc_shimalik"]
    supported_methods_embeddings = ["unnorm_ssc_hope", "norm_ssc_hope", "unnorm_hsc_hope", "norm_hsc_hope", "hsc_shimalik_hope", "ssc_shimalik_hope",
                                       "unnorm_ssc_laplacianeigenmaps", "norm_ssc_laplacianeigenmaps", "unnorm_hsc_laplacianeigenmaps", "norm_hsc_laplacianeigenmaps",
                                       "hsc_shimalik_laplacianeigenmaps", "ssc_shimalik_laplacianeigenmaps"]
    supported_methods_communities = ["louvain", "leiden", "walktrap", "asyn_fluidc", "greedy_modularity", "infomap", "spinglass"]
    
    if platform == "win32":
        # Windows...
        supported_methods_communities.remove("infomap")

    if ((method in supported_methods_clustering)|(method in supported_methods_communities)|(method in supported_methods_embeddings)):
    
        if method in supported_methods_clustering:
            name = "Clusters"
        if method in supported_methods_communities:
            name = "Communities"
        if method in supported_methods_embeddings:
            name = "ClustersEmbeddings"

    residue_names_0 = np.array(residue_names[:, 0], dtype = int)
    residue_names_1 = np.array(residue_names[:, 1], dtype = str)

    labels_u, counts = np.unique(labels, return_counts=True)

    for i in range(len(labels_u)):
    
        print(labels_u[i], counts[i])
       
    dict_node_cluster_0 = dict ()   
    for i, label in enumerate(labels):
    
      dict_node_cluster_0[str (residue_names_0[i])] = label    
    
    residue_clusters = np.array(list (dict_node_cluster_0.items()))
    residue_names_cluster = dict ()   
    
    k = int (max(labels)) + 1
    
    for label in range(k):
    
        temp = []
        
        for (resi, res_name) in residue_names:

            if ((residue_names[int (resi), 0]==residue_clusters[int (resi), 0]) and (str (int ((residue_clusters[int (resi), 1]))) == str (label))): 
            
                temp.append(res_name)
        
        print(len(temp))
        residue_names_cluster[label] = temp
        print(f"{name} {label}: ", residue_names_cluster[label]) 
    
    #TO DELETE(maybe)
    
    summary_output_path = "{}{}{}{}".format(output_path, method, add_slash_to_path, "Summary")
    
    if (not os.path.exists(summary_output_path)):
    
        os.makedirs(summary_output_path)
        
    f = open("{}{}{}_{}_{}_k{}.txt".format(summary_output_path, add_slash_to_path, p_name, "{}_Summary".format(name), method, k),"w")
    
    for label in range(k):
    
        f.write("{} {}: {} ".format(name, label, residue_names_cluster[label]))
        f.write(" ")
    
    f.close()
    #END TO DELETE 
            
    dict_node_cluster_1 = dict ()
    
    for i, label in enumerate(labels):
    
        dict_node_cluster_1[residue_names_1[i]] = int (label)    
        
    method_output_path = "{}{}{}{}".format(output_path, method, add_slash_to_path, name)

    if (not os.path.exists(method_output_path)):
    
        os.makedirs(method_output_path)
        
    if (name == "Clusters"):
    
        f = open("{}{}{}_{}_{}_k{}.txt".format(method_output_path, add_slash_to_path, p_name, name, method, k),"w")
        f.write(str (dict_node_cluster_1))
        f.close()
        return dict_node_cluster_1
        
    elif (name == "ClustersEmbeddings"):
        
        if beta is not None:

            f = open("{}{}{}_{}_{}_d{}_beta{}_k{}.txt".format(method_output_path, add_slash_to_path, p_name, name, method, d, beta, k),"w")

        else:
            
            f = open("{}{}{}_{}_{}_d{}_k{}.txt".format(method_output_path, add_slash_to_path, p_name, name, method, d, k),"w")
            
        f.write(str (dict_node_cluster_1))
        f.close()
        return dict_node_cluster_1
        
    elif (name == "Communities"):
    
        f = open("{}{}{}_{}_{}_ncoms{}.txt".format(method_output_path, add_slash_to_path, p_name, name, method, k),"w")    
        f.write(str (dict_node_cluster_1))
        f.close()
        
        return dict_node_cluster_1
        
    else:
        raise Exception ("method {} not supported".format(method))
        

#COMMUNITY EXTRACTION

def extract_labels_from_coms(num_nodes,coms, algorithm_name):
  
    if (algorithm_name != "Asyn FluidC"):
        n_coms = len(coms) #numero di comunità 
        print("number of {} communities: {}".format(algorithm_name, n_coms))
    
    labels = np.zeros((num_nodes, 1))
    dict_node_algorithm_com = dict ()

    for label, com in enumerate(coms):
  
        for i, node in enumerate(com):
    
            node =  int (float (node))
            dict_node_algorithm_com[node] = label
            labels[node] = label
  
    return labels

def louvain(G):

    louvain = algorithms.louvain(G)
    coms = louvain.communities
    num_nodes = G.number_of_nodes()
    labels = extract_labels_from_coms(num_nodes, coms, "Louvain")

    return labels

def leiden(G):

    leiden = algorithms.leiden(G)
    coms = leiden.communities
    num_nodes = G.number_of_nodes()
    labels = extract_labels_from_coms(num_nodes, coms, "Leiden")

    return labels  
 
def walktrap(G):

    walktrap = algorithms.walktrap(G)
    coms = walktrap.communities
    num_nodes = G.number_of_nodes()
    labels = extract_labels_from_coms(num_nodes, coms, "Walktrap")
      
    return labels
  
def greedy_modularity(G):

    coms = greedy_modularity_communities(G)
    num_nodes = G.number_of_nodes()
    labels = extract_labels_from_coms(num_nodes, coms, "Greedy Modularity")

    return labels

def infomap(G):
    
    coms = algorithms.infomap(G).communities
    num_nodes = G.number_of_nodes()
    labels = extract_labels_from_coms(num_nodes, coms, "Infomap")

    return labels

def asyn_fluidc(G, k):
  
    coms = asyn_fluidc_(G, k)
    num_nodes = G.number_of_nodes()
    print("number of Asyn FluidC communities: {}".format(k))
    labels = extract_labels_from_coms(num_nodes, coms, "Asyn FluidC")
  
    return labels
  
def spinglass(G):

    coms = algorithms.spinglass(G).communities
    num_nodes = G.number_of_nodes()
    labels = extract_labels_from_coms(num_nodes, coms, "Spin Glass")
  
    return labels
  
#SPECTRAL CLUSTERING

def degree_matrix(A):
  
    n = A.shape[0]
    diag = np.zeros((n, n)) 
    rows, cols = A.nonzero()
    
    for row, col in zip(rows, cols):
    
        diag[row, row] += 1

    return diag

def compute_laplacian_matrix(A):

    D = degree_matrix(A)
    L = D-A
    return L
 
def compute_normalized_laplacian(A):

    D = degree_matrix(A)
    D_inv_sqrt = np.linalg.inv(np.sqrt(D))
    L_sym = 1 - np.dot(D_inv_sqrt, A).dot(D_inv_sqrt)
    return L_sym

def computeSortEigens(mat):

    eigenvalues, eigenvectors = scipy.linalg.eig(mat)
    
    idx = eigenvalues.argsort()
    sortedEigenvalues = eigenvalues[idx].real
    sortedEigenvectors = eigenvectors[:,idx].real
  
    return sortedEigenvalues, sortedEigenvectors
  
def computeBestK(eigenvalues, n_k = 1):

    sortedEigenvalues = sorted(eigenvalues)[::-1]
    max_eigengap = 0              
    n = len(sortedEigenvalues)
    all_k = np.arange(0, n)
    best_k = 0
    all_eigengaps = np.zeros((n))
  
    for index, k in enumerate(all_k):
    
        if (k!=0):

            eigengap = abs(sortedEigenvalues[k]-sortedEigenvalues[k-1])
            all_eigengaps[k] = eigengap

            if (eigengap > max_eigengap):
        
                best_k = k
                max_eigengap = eigengap

        else: continue
  
    idx = all_eigengaps.argsort()[::-1]
    real_best_ks = all_k[idx]
  
    if ((best_k > 60)|(best_k==1)):
  
        best_k = 0
        for k in real_best_ks:
    
            if ((k<60)&(k>1)):
      
                best_k = k
                break

    best_ks = real_best_ks[(real_best_ks<60) & (real_best_ks>1)][:n_k]
    print("Best k: ", best_k)
    #print("Real Best {} ks: ".format(n_k), real_best_ks)
    print("Choosen Best {} ks: ".format(n_k), best_ks)
  
    return best_ks

def hardSpectralClustering(A, n_clusters = None, norm=False, embedding=None, d=2, beta=0.01):
    
    supported_embeddings = ["HOPE", "LaplacianEigenmaps"]
    
    if embedding is None:
    
        if norm:
            L = compute_normalized_laplacian(A)
        else:
            L = compute_laplacian_matrix(A)

        sortedEigenvalues, sortedEigenvectors = computeSortEigens(L)             #sorted eigenvalues/eigenvectors
        train = sortedEigenvectors[:, :n_clusters]
                   
    else:
    
        if (embedding in supported_embeddings):
    
            if (embedding == "HOPE"): #embedding dimension (d) and decay factor (beta) as inputs
                model = HOPE(d=d, beta=beta)
            if (embedding == "node2vec"):
                #model = node2vec(d=d, max_iter=1, walk_len=80, num_walks=10, con_size=10, ret_p=1, inout_p=1)
                pass
            if (embedding == "LaplacianEigenmaps"):
                model = LaplacianEigenmaps(d=d)
            train, t = model.learn_embedding(graph=nx.from_numpy_matrix(A))
        else: raise Exception ("embedding {} not supported".format(embedding))
    
    km = KMeans(n_clusters=n_clusters).fit(train)
    labels = km.labels_
 
    return labels

def softSpectralClustering(A, n_clusters = None, norm=False, embedding = None,  d=2, beta=0.01):
     
    supported_embeddings = ["HOPE", "LaplacianEigenmaps"]
    
    if embedding is None:   
        
        if norm:
            L = compute_normalized_laplacian(A)
        else:
            L = compute_laplacian_matrix(A)
            
        sortedEigenvalues, sortedEigenvectors = computeSortEigens(L)            #sorted eigenvalues/eigenvectors
        train = sortedEigenvectors[:, :n_clusters]
                
    else:
    
        if (embedding in supported_embeddings):
             
            if (embedding == "HOPE"):
                model = HOPE(d=d, beta=beta)
            if (embedding == "node2vec"):
                #model = node2vec(d=d, max_iter=1, walk_len=80, num_walks=10, con_size=10, ret_p=1, inout_p=1)
                pass
            if (embedding == "LaplacianEigenmaps"):
                model = LaplacianEigenmaps(d=d)
            train, t = model.learn_embedding(graph=nx.from_numpy_matrix(A))
        else: raise Exception ("embedding {} not supported".format(embedding))
    
    fcm = FCM(n_clusters=n_clusters)
    fcm.fit(train)
    labels = fcm.predict(train)

    return labels
  
def ssc_shimalik(A, n_clusters = None, embedding=None, d=2, beta=0.01):
    
    supported_embeddings = ["HOPE", "LaplacianEigenmaps"]
    
    if embedding is None:
        
        L = compute_laplacian_matrix(A)
        D = degree_matrix(A)
        eigenvalues, eigenvectors = eigh(L, D, eigvals_only=False)
        idx = eigenvalues.argsort()
        sortedEigenvalues = eigenvalues[idx].real
        sortedEigenvectors = eigenvectors[:,idx].real
        train = sortedEigenvectors[:, :n_clusters]
       
    else:
    
        if (embedding in supported_embeddings):
    
            if (embedding == "HOPE"):
                model = HOPE(d=d, beta=beta)
            if (embedding == "node2vec"):
                #model = node2vec(d=d, max_iter=1, walk_len=80, num_walks=10, con_size=10, ret_p=1, inout_p=1)
                pass
            if (embedding == "LaplacianEigenmaps"):
                model = LaplacianEigenmaps(d=d)
            train, t = model.learn_embedding(graph=nx.from_numpy_matrix(A))
        else: raise Exception ("embedding {} not supported".format(embedding))
        
    fcm = FCM(n_clusters=n_clusters)
    fcm.fit(train)
    labels = fcm.predict(train)

    return labels

def hsc_shimalik(A, n_clusters = None, embedding = None, d=2, beta=0.01):

    supported_embeddings = ["HOPE", "LaplacianEigenmaps"]
  
    if embedding is None:
  
        L = compute_laplacian_matrix(A)
        D = degree_matrix(A)
        eigenvalues, eigenvectors = eigh(L, D, eigvals_only=False)
        idx = eigenvalues.argsort()
        sortedEigenvalues = eigenvalues[idx].real
        sortedEigenvectors = eigenvectors[:,idx].real
        
        if n_clusters is None:
            n_clusters = computeBestK(sortedEigenvalues)
        
        train = sortedEigenvectors[:, :n_clusters]
        
    else:
  
        if (embedding in supported_embeddings):
            if (embedding == "HOPE"):
              model = HOPE(d=d, beta=beta)
            if (embedding == "node2vec"):
              #model = node2vec(d=2, max_iter=1, walk_len=80, num_walks=10, con_size=10, ret_p=1, inout_p=1)
              pass
            if (embedding == "LaplacianEigenmaps"):
              model = LaplacianEigenmaps(d=d)
            train, t = model.learn_embedding(graph=nx.from_numpy_matrix(A))
        else:
            raise Exception ("embedding {} not supported".format(embedding))

    km = KMeans(n_clusters=n_clusters).fit(train)
    labels = km.labels_

    return labels
  
 #spectral clustering 
def unnorm_hsc(A, n_clusters = None, norm=False, embedding=None, d=None, beta=None):

    labels = hardSpectralClustering(A, n_clusters, norm, embedding, d, beta)
    return labels
    
def norm_hsc(A, n_clusters = None, norm=True, embedding=None, d=None, beta=None):

    labels = hardSpectralClustering(A, n_clusters, norm, embedding, d, beta)
    return labels
    
def unnorm_ssc(A, n_clusters = None, norm=False, embedding=None, d=None, beta=None):

    labels = softSpectralClustering(A, n_clusters, norm, embedding, d, beta)
    return labels

def norm_ssc(A, n_clusters = None, norm=True, embedding=None, d=None, beta=None):

    labels = softSpectralClustering(A, n_clusters, norm, embedding, d, beta)
    return labels

#EMBEDDINGS
def unnorm_hsc_hope(A, n_clusters = None, norm=False, embedding="HOPE", d=2, beta=0.01):

    labels = hardSpectralClustering(A, n_clusters, norm, embedding, d, beta)
    return labels
    
def norm_hsc_hope(A, n_clusters = None, norm=True, embedding="HOPE", d=2, beta=0.01):

    labels = hardSpectralClustering(A, n_clusters, norm, embedding, d, beta)
    return labels

def unnorm_ssc_hope(A, n_clusters = None, norm=False, embedding="HOPE", d=2, beta=0.01):

    labels = softSpectralClustering(A, n_clusters, norm, embedding, d, beta)
    return labels
    
def norm_ssc_hope(A, n_clusters = None, norm=True, embedding="HOPE", d=2, beta=0.01):

    labels = softSpectralClustering(A, n_clusters, norm, embedding, d, beta)
    return labels
    
def hsc_shimalik_hope(A, n_clusters = None, embedding="HOPE", d=2, beta=0.01):

    labels = hsc_shimalik(A, n_clusters, embedding, d, beta)
    return labels
    
def ssc_shimalik_hope(A, n_clusters = None, embedding="HOPE", d=2, beta=0.01):

    labels = ssc_shimalik(A, n_clusters, embedding, d, beta)
    return labels

def unnorm_hsc_laplacianeigenmaps(A, n_clusters = None, norm=False, embedding="LaplacianEigenmaps", d=2, beta=None):

    labels = hardSpectralClustering(A, n_clusters, norm, embedding, d, beta)
    return labels
    
def norm_hsc_laplacianeigenmaps(A, n_clusters = None, norm=True, embedding="LaplacianEigenmaps", d=2, beta=None):

    labels = hardSpectralClustering(A, n_clusters, norm, embedding, d, beta)
    return labels
    
def unnorm_ssc_laplacianeigenmaps(A, n_clusters = None, norm=False, embedding="LaplacianEigenmaps", d=2, beta=None):

    labels = softSpectralClustering(A, n_clusters, norm, embedding, d, beta)
    return labels
    
def norm_ssc_laplacianeigenmaps(A, n_clusters = None, norm=True, embedding="LaplacianEigenmaps", d=2, beta=None):


    labels = softSpectralClustering(A, n_clusters, norm, embedding, d, beta)
    return labels
    
def hsc_shimalik_laplacianeigenmaps(A, n_clusters = None, embedding="LaplacianEigenmaps", d=2, beta=None):

    labels = hsc_shimalik(A, n_clusters, embedding, d, beta)
    return labels
    
def ssc_shimalik_laplacianeigenmaps(A, n_clusters = None, embedding="LaplacianEigenmaps", d=2, beta=None):

    labels = ssc_shimalik(A, n_clusters, embedding, d, beta)
    return labels
    
#NEW

def betweenness(G, p, residue_names, n):

    bc = betweenness_centrality(G)
    bc = {int (float (k)):v for k,v in bc.items()}
    residue_names_1 = np.array(residue_names[:, 1], dtype = str)
    dict_node_centrality = dict ()
    
    for i, cent in enumerate(bc.values()):
    
        dict_node_centrality[residue_names_1[i]] = cent    
        
    sorted_bc = sorted(dict_node_centrality.items(), key=itemgetter(1), reverse=True)
    print("Top {} nodes by betweenness centrality".format(n))
    for d in sorted_bc[:n]:
        print(d)
    
    return dict_node_centrality

def eigenvector(G, p, residue_names, n):

    ec = eigenvector_centrality(G, max_iter=10000)
    ec = {int (float (k)):v for k,v in ec.items()}
    residue_names_1 = np.array(residue_names[:, 1], dtype = str)
    dict_node_centrality = dict ()
    
    for i, cent in enumerate(ec.values()):
    
        dict_node_centrality[residue_names_1[i]] = cent      
        
    sorted_ec = sorted(dict_node_centrality.items(), key=itemgetter(1), reverse=True)
    print("Top {} nodes by eigenvector centrality".format(n))
    for d in sorted_ec[:n]:
        print(d)
        
    return dict_node_centrality

def degree_c(G, p, residue_names, n):

    dc = degree_centrality(G)
    dc = {int (float (k)):v for k,v in dc.items()}
    residue_names_1 = np.array(residue_names[:, 1], dtype = str)
    dict_node_centrality = dict ()
    
    for i, cent in enumerate(dc.values()):
    
        dict_node_centrality[residue_names_1[i]] = cent     
        
    sorted_dc = sorted(dict_node_centrality.items(), key=itemgetter(1), reverse=True)
    print("Top {} nodes by degree centrality".format(n))
    for d in sorted_dc[:n]:
        print(d)
    
    return dict_node_centrality

def closeness(G, p, residue_names, n):

    cc = closeness_centrality(G)
    cc = {int (float (k)):v for k,v in cc.items()}
    residue_names_1 = np.array(residue_names[:, 1], dtype = str)
    dict_node_centrality = dict ()
    
    for i, cent in enumerate(cc.values()):
    
        dict_node_centrality[residue_names_1[i]] = cent  
        
    sorted_cc = sorted(dict_node_centrality.items(), key=itemgetter(1), reverse=True)
    print("Top {} nodes by closeness_centrality".format(n))
    for d in sorted_cc[:n]:
        print(d)
        
    return dict_node_centrality
    
def partecipation_coefs(G, clusters):
  
    A = nx.to_numpy_matrix(G)
    n = A.shape[0]
    P = dict()
    k_s = np.zeros((n))

    for i in range(n):
        k_i = np.sum(A[i,:]) 
        k_si = 0
    
        for j in range(n):
            if (i!=j):
                if ((clusters[i] == clusters[j]) and (A[i,j]!=0)):#se il nodo i e il nodo j sono dello stesso cluster e c'è un arco che li connette
                    k_si += A[i,j]

        k_s[i] = k_si
        P[i] = 1 - (k_s[i]/k_i)**2

    return P

def z_score(G, clusters):

    A = nx.to_numpy_matrix(G)
    n = A.shape[0]
    z = dict()
    k_s = np.zeros((n))

    for i in range(n):
        k_i = np.sum(A[i,:]) 
        k_si = 0
        
        for j in range(n):
            if (i!=j):
                if ((clusters[i] == clusters[j]) and (A[i,j]!=0)):#se il nodo i e il nodo j sono dello stesso cluster e c'è un arco che li connette
                    k_si += A[i,j]

        k_s[i] = k_si
      
        mean_k_si = np.mean(k_s)
        sigma_k_si = np.std(k_s)
        for i in range(n): 
            k_i = np.sum(A[i,:]) 
            z[i] = (k_i - mean_k_si) / sigma_k_si
  
    return z

def plot_z_p(G, clusters):
    
    part_coefs = partecipation_coefs(G, clusters)
    z_scores = z_score(G, clusters)
    plt.scatter(np.array(list(part_coefs.values())).astype(float), np.array(list(z_scores.values())).astype(float), facecolors='none', edgecolors='b')
    plt.xlabel("part_coefs")
    plt.ylabel("z_scores")
    
def color_map_clustering(clusters):
    
    n = len(clusters)
    labels_unique, counts = np.unique(clusters, return_counts=True)
        
    idx = labels_unique.argsort()
    labels_unique = labels_unique[idx]
    counts = counts[idx]
    
    cluster_map = np.zeros((n,n))
    
    j_to_start = 0
    for idx, label_unique in enumerate(labels_unique): 
        
        temp_count = 0
        count_label = counts[idx]  
                 
        while(temp_count < count_label):
            for i, label in enumerate(clusters):
    
                if label == label_unique:
                
                    if(j_to_start+count_label<n):
                        cluster_map[i][j_to_start:j_to_start+count_label] = label+1
                    else:
                        cluster_map[i][j_to_start:] = label+1
                    temp_count += 1

                if temp_count >= count_label: 
                    j_to_start += temp_count
                    break                 
                
    for i in range(n):
        cluster_map[i][i]=0
                
    print(cluster_map)
    plt.figure(figsize=(8,8))
    plt.matshow(cluster_map, vmin=0, vmax=max(clusters)+1, cmap='cividis', fignum=1)
    plt.legend(labels=labels_unique)
    plt.show()
