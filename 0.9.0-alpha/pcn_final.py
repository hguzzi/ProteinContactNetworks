#IMPORT LIBRARIES
from prody import *
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

from gem.embedding.hope     import HOPE
from gem.embedding.lap      import LaplacianEigenmaps
from gem.embedding.node2vec import node2vec

from cdlib import algorithms
from networkx.algorithms.community.centrality import girvan_newman
from networkx.algorithms.community.asyn_fluid import asyn_fluidc
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
    
def checkIfFilesExists(files, initial_choice, directory):
    
    not_existing_files = []
    all_files_exists = True
    
    for file in files:
        file_path = "{}{}".format(directory, file)
        if((not os.path.isfile(file_path))):
            all_files_exists = False
            not_existing_files.append(file)
    
    if(not all_files_exists):
        if(initial_choice == ""):
            for protein in not_existing_files:
                urllib.request.urlretrieve('http://files.rcsb.org/download/{}.pdb'.format(protein), '{}{}.pdb'.format(directory, protein))
        else:
            raise Exception("Files {} not exists in directory '{}'.".format(str(not_existing_files), directory))
        
    return all_files_exists    
    
def readPDBFile(pbdFilePath):
  atoms = []
  with open(pbdFilePath) as pdbfile:
      for line in pdbfile:
          if line[:4] == 'ATOM':
              # Split the line
                splitted_line = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54], line[56:61], line[62:66]]
                atoms.append(splitted_line)
                
                # To format again the pdb file with the fields extracted
                #print("%-6s%5s %4s %3s %s %4s    %8s%8s%8s   %3s%3s\n"%tuple(splitted_line))
  return np.array(atoms)

def getResidueDistance(atoms):
  
  residues = []
  residues_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 
                   'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 
                   'TRP', 'TYR']
  
  dropped = []

  last_residue_num = 0
  last_chain = ''

  for i, atom in enumerate(atoms):

    residue_name = atom[3]
    residue_chain = atom[4]
    residue_num = int (atom[5])
    
    if residue_name in residues_list:   

      if (residue_num != last_residue_num):

        if (atom[2]==' C  '):
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
 
def getResiduesSequence(pbdFilePath):
  
  seq_res = []
  residues_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
  
  with open(pbdFilePath) as pdbfile:
    for line in pdbfile:
      if line[:6]=='SEQRES':
        splitted_line = [line[19:22], line[23:26], line[27:30], line[31:34], line[35:38],
                        line[39:42], line[43:46], line[47:50], line[51:54], line[55:58],
                        line[59:62], line[63:66], line[67:70]]
        for residue in splitted_line:
          if residue in residues_list:
            seq_res.append(residue)
  return np.array(seq_res)

def read_adj_mat(adj_filepath, p):
    if (os.path.isfile("{}{}_adj_mat.txt".format(adj_filepath, p))):
        adj = np.loadtxt("{}{}_adj_mat.txt".format(adj_filepath, p))
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
          #print(i+'out of '+n)
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
    np.savetxt("{}Edgelists{}{}_edgelist.txt".format(output_path, add_slash_to_path, p), np.array(edge_list), fmt='%.2f')
    print("saved edge list")
        
    if not os.path.exists("{}Adj".format(output_path)):
        os.makedirs("{}Adj".format(output_path))
    np.savetxt("{}Adj{}{}_adj_mat.txt".format(output_path, add_slash_to_path, p), adj, fmt='%.2f')
    print("saved adj matrix")
    print("adj shape: ", adj.shape)
    
    return adj

#SAVE LABELS
def save_labels(output_path, labels, residue_names, p_name, method=None, d=None, beta=None):

    supported_methods_clustering = ["Unnorm_SSC", "Norm_SSC", "Unnorm_HSC", "Norm_HSC", "SKL_Spectral", "HSC_ShiMalik", "SSC_ShiMalik"]
    supported_methods_embeddings = ["Unnorm_SSC_HOPE", "Norm_SSC_HOPE", "Unnorm_HSC_HOPE", "Norm_HSC_HOPE", "HSC_ShiMalik_HOPE", "SSC_ShiMalik_HOPE",
                                       "Unnorm_SSC_LaplacianEigenmaps", "Norm_SSC_LaplacianEigenmaps", "Unnorm_HSC_LaplacianEigenmaps", "Norm_HSC_LaplacianEigenmaps",
                                       "HSC_ShiMalik_LaplacianEigenmaps", "SSC_ShiMalik_LaplacianEigenmaps"]
    supported_methods_communities = ["Louvain", "Leiden", "Walktrap", "Asyn_fluidc", "Greedy_modularity", "Infomap", "Spinglass"]

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
            
    dict_node_cluster_1 = dict ()
    for i, label in enumerate(labels):
        dict_node_cluster_1[residue_names_1[i]] = int (label)    

    if (not os.path.exists("{}{}{}{}".format(output_path, method, add_slash_to_path, name))):
        os.makedirs("{}{}{}{}".format(output_path, method, add_slash_to_path, name))
        
    if (name == "Clusters"):
        f = open("{}{}{}{}{}{}_{}_{}_k{}.txt".format(output_path, method, add_slash_to_path, name, add_slash_to_path, p_name, name, method, k),"w")
        f.write(str (dict_node_cluster_1))
        f.close()
        return dict_node_cluster_1
        
    elif (name == "ClustersEmbeddings"):
        
        if beta is not None:

            f = open("{}{}{}{}{}{}_{}_{}_d{}_beta{}_k{}.txt".format(output_path, method, add_slash_to_path, name, add_slash_to_path, p_name, name, method, d, beta, k),"w")

        else:
            
            f = open("{}{}{}{}{}{}_{}_{}_d{}_k{}.txt".format(output_path, method, add_slash_to_path, name, add_slash_to_path, p_name, name, method, d, k),"w")
            
        f.write(str (dict_node_cluster_1))
        f.close()
        return dict_node_cluster_1
        
    elif (name == "Communities"):
        f = open("{}{}{}{}{}{}_{}_{}_k{}.txt".format(output_path, method, add_slash_to_path, name, add_slash_to_path, p_name, name, method, k),"w")    
        f.write(str (dict_node_cluster_1))
        f.close()
        
        return dict_node_cluster_1
        
    else:
        raise Exception ("method {} not supported".format(method))
#COMMUNITY EXTRACTION
def Louvain(G):
  print('Start Computation of Louvain Algorithm')
  louvain = algorithms.louvain(G)
  coms = louvain.communities
  n_coms = len(coms) #numero di comunità 
  print("number of Louvain communities: {}".format(n_coms))
  
  labels = np.zeros((G.number_of_nodes(), 1))
  dict_node_Louvain_com = dict ()

  for label, com in enumerate(coms):
    for i, node in enumerate(com):
      node =  int (float (node))
      dict_node_Louvain_com[node] = label
      labels[node] = label
  
  labels_u, counts = np.unique(labels, return_counts=True)

  return labels

def Leiden(G):

  leiden = algorithms.leiden(G)
  coms = leiden.communities
  n_coms = len(coms) #numero di comunità 
  print("number of Leiden communities: {}".format(n_coms))
  
  labels = np.zeros((G.number_of_nodes(), 1))
  dict_node_Leiden_com = dict ()

  for label, com in enumerate(coms):
    for i, node in enumerate(com):
      node =  int (float (node))
      dict_node_Leiden_com[node] = label
      labels[node] = label

  labels_u, counts = np.unique(labels, return_counts=True)
  return labels  
 
def Walktrap(G):

  walktrap = algorithms.walktrap(G)
  coms = walktrap.communities
  n_coms = len(coms) #numero di comunità 
  print("number of Walktrap communities: {}".format(n_coms))
  
  labels = np.zeros((G.number_of_nodes(), 1))
  dict_node_Walktrap_com = dict ()

  for label, com in enumerate(coms):
    for i, node in enumerate(com):
      node =  int (float (node))
      dict_node_Walktrap_com[node] = label
      labels[node] = label
      
  labels_u, counts = np.unique(labels, return_counts=True)

  return labels
  
def Greedy_modularity(G):
  coms = greedy_modularity_communities(G)
  dict_node_greedy_modularity_com = dict ()
  labels = np.zeros((G.number_of_nodes(), 1))
  
  for label, com in enumerate(coms):
    for i, node in enumerate(com):
      node =  int (float (node))
      dict_node_greedy_modularity_com[node] = label
      labels[node] = label

  labels_u, counts = np.unique(labels, return_counts=True)
  print("number of greedy modularity communities: {}".format(int (max(labels_u)+1)))
  return labels

def Infomap(G):
    
    coms = algorithms.infomap(G).communities
        
    dict_node_infomap_com = dict ()
    labels = np.zeros((G.number_of_nodes(), 1))
  
    for label, com in enumerate(coms):
        for i, node in enumerate(com):
            node =  int (float (node))
            dict_node_infomap_com[node] = label
            labels[node] = label

    labels_u, counts = np.unique(labels, return_counts=True)
    print("number of infomap communities: {}".format(int (max(labels_u)+1)))
    return labels

def Asyn_fluidc(G):
  
  k = int( input("Insert k for Asyn_fluidc algorithm: "))
  print("number of Asyn_fluidc communities: ", k)
  coms = asyn_fluidc(G, k)

  dict_node_Asyn_fluidc_com = dict ()
  labels = np.zeros((G.number_of_nodes(), 1))
  
  for label, com in enumerate(coms):
    for i, node in enumerate(com):
      node =  int (float (node))
      dict_node_Asyn_fluidc_com[node] = label
      labels[node] = label

  labels_u, counts = np.unique(labels, return_counts=True)
  return labels
  
def Spinglass(G):
  coms = algorithms.spinglass(G).communities
  print(str (coms))
  dict_node_spinglass_com = dict ()
  labels = np.zeros((G.number_of_nodes(), 1))
  
  for label, com in enumerate(coms):
    for i, node in enumerate(com):
      node =  int (float (node))
      dict_node_spinglass_com[node] = label
      labels[node] = label

  labels_u, counts = np.unique(labels, return_counts=True)
  print("number of spinglass communities: {}".format(int (max(labels_u)+1)))
  return labels
  
#SPECTRAL CLUSTERING

def degree_matrix(A):
  
    n = A.shape[0]
    diag = np.zeros((n, n)) # basically dimensions of your graph
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
  
def computeBestK(eigenvalues, n_k = 2):

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
  
  if ((best_k > 27)|(best_k==1)):
    best_k = 0
    for k in real_best_ks:
      if ((k<27)&(k>1)):
        best_k = k
        break

  best_ks = real_best_ks[(real_best_ks<27) & (real_best_ks>1)][:n_k]
  print("Best k: ", best_k)
  print("Real Best {} ks: ".format(n_k), real_best_ks)
  print("Choosen Best {} ks: ".format(n_k), best_ks)
  
  return best_ks

def hardSpectralClustering(A, n_clusters = None, plot_eig=False, plot_fiedler=False, norm=False, embedding=None, d=2, beta=0.01):
    
    supported_embeddings = ["node2vec", "HOPE", "LaplacianEigenmaps"]
    if embedding is None:
        if norm:
            L = compute_normalized_laplacian(A)
        else:
            L = compute_laplacian_matrix(A)

        sortedEigenvalues, sortedEigenvectors = computeSortEigens(L)             #sorted eigenvalues/eigenvectors

        if n_clusters is None:
            n_clusters = computeBestK(sortedEigenvalues, n_k = 5)
        
        j, fiedler_eigenvalue, fiedler_eigenvector = next(([j, eigenvalue, sortedEigenvectors[j]] for j, eigenvalue in enumerate(sortedEigenvalues) if (round(eigenvalue, 6) != 0)), None)
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

def softSpectralClustering(A, n_clusters = None, plot_eig=False, plot_fiedler=False, norm=False, embedding = None,  d=2, beta=0.01):
     
    supported_embeddings = ["node2vec", "HOPE", "LaplacianEigenmaps"]
    if embedding is None:   
        
        if norm:
            L = compute_normalized_laplacian(A)
        else:
            L = compute_laplacian_matrix(A)

        sortedEigenvalues, sortedEigenvectors = computeSortEigens(L)            #sorted eigenvalues/eigenvectors
        j, fiedler_eigenvalue, fiedler_eigenvector = next(([j, eigenvalue, sortedEigenvectors[j]] for j, eigenvalue in enumerate(sortedEigenvalues) if (round(eigenvalue, 6) != 0)), None)
        train = sortedEigenvectors[:, :n_clusters]
        
        if n_clusters is None:
            n_clusters = computeBestK(sortedEigenvalues,  n_k = 5)
            
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
    centers = fcm.centers
    print("labels: ", labels)
    print("centers: ", centers)

    return labels
  
def SSC_ShiMalik(A, n_clusters = None, plot_eig=False, plot_fiedler=False, embedding=None, d=2, beta=0.01):
    
    supported_embeddings = ["node2vec", "HOPE", "LaplacianEigenmaps"]
    if embedding is None:
        L = compute_laplacian_matrix(A)
        D = degree_matrix(A)
        eigenvalues, eigenvectors = eigh(L, D, eigvals_only=False)
        idx = eigenvalues.argsort()
        sortedEigenvalues = eigenvalues[idx].real
        sortedEigenvectors = eigenvectors[:,idx].real
        j, fiedler_eigenvalue, fiedler_eigenvector = next(([j, eigenvalue, sortedEigenvectors[j]] for j, eigenvalue in enumerate(sortedEigenvalues) if (round(eigenvalue, 6) != 0)), None)

        if n_clusters is None:
            n_clusters = computeBestK(sortedEigenvalues)

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
    centers = fcm.centers
    print("labels: ", labels)
    print("centers: ", centers)

    return labels

def HSC_ShiMalik(A, n_clusters = None, plot_eig=False, plot_fiedler=False, embedding = None, d=2, beta=0.01):

  supported_embeddings = ["node2vec", "HOPE", "LaplacianEigenmaps"]
  if embedding is None:
    L = compute_laplacian_matrix(A)
    D = degree_matrix(A)
    eigenvalues, eigenvectors = eigh(L, D, eigvals_only=False)
    idx = eigenvalues.argsort()
    sortedEigenvalues = eigenvalues[idx].real
    sortedEigenvectors = eigenvectors[:,idx].real
    j, fiedler_eigenvalue, fiedler_eigenvector = next(([j, eigenvalue, sortedEigenvectors[j]] for j, eigenvalue in enumerate(sortedEigenvalues) if (round(eigenvalue, 6) != 0)), None)
    train = sortedEigenvectors[:, :n_clusters]
    
    if n_clusters is None:
      n_clusters = computeBestK(sortedEigenvalues)
  
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
  
 #
def Unnorm_HSC(A, n_clusters = None, plot_eig=False, plot_fiedler=False, norm=False, embedding=None, d=None, beta=None):
    labels = hardSpectralClustering(A, n_clusters, plot_eig, plot_fiedler, norm, embedding, d, beta)
    return labels
    
def Norm_HSC(A, n_clusters = None, plot_eig=False, plot_fiedler=False, norm=True, embedding=None, d=None, beta=None):
    labels = hardSpectralClustering(A, n_clusters, plot_eig, plot_fiedler, norm, embedding, d, beta)
    return labels
    
def Unnorm_SSC(A, n_clusters = None, plot_eig=False, plot_fiedler=False, norm=False, embedding=None, d=None, beta=None):
    labels = softSpectralClustering(A, n_clusters, plot_eig, plot_fiedler, norm, embedding, d, beta)
    return labels

def Norm_SSC(A, n_clusters = None, plot_eig=False, plot_fiedler=False, norm=True, embedding=None, d=None, beta=None):
    labels = softSpectralClustering(A, n_clusters, plot_eig, plot_fiedler, norm, embedding, d, beta)
    return labels

#EMBEDDINGS
def Unnorm_HSC_HOPE(A, n_clusters = None, plot_eig=False, plot_fiedler=False, norm=False, embedding="HOPE", d=2, beta=0.01):
    labels = hardSpectralClustering(A, n_clusters, plot_eig, plot_fiedler, norm, embedding, d, beta)
    return labels
    
def Norm_HSC_HOPE(A, n_clusters = None, plot_eig=False, plot_fiedler=False, norm=True, embedding="HOPE", d=2, beta=0.01):
    labels = hardSpectralClustering(A, n_clusters, plot_eig, plot_fiedler, norm, embedding, d, beta)
    return labels

def Unnorm_SSC_HOPE(A, n_clusters = None, plot_eig=False, plot_fiedler=False, norm=False, embedding="HOPE", d=2, beta=0.01):
    labels = softSpectralClustering(A, n_clusters, plot_eig, plot_fiedler, norm, embedding, d, beta)
    return labels
    
def Norm_SSC_HOPE(A, n_clusters = None, plot_eig=False, plot_fiedler=False, norm=True, embedding="HOPE", d=2, beta=0.01):
    labels = softSpectralClustering(A, n_clusters, plot_eig, plot_fiedler, norm, embedding, d, beta)
    return labels
    
def HSC_ShiMalik_HOPE(A, n_clusters = None, plot_eig=False, plot_fiedler=False, embedding="HOPE", d=2, beta=0.01):
    labels = HSC_ShiMalik(A, n_clusters, plot_eig, plot_fiedler, embedding, d, beta)
    return labels
    
def SSC_ShiMalik_HOPE(A, n_clusters = None, plot_eig=False, plot_fiedler=False, embedding="HOPE", d=2, beta=0.01):
    labels = SSC_ShiMalik(A, n_clusters, plot_eig, plot_fiedler, embedding, d, beta)
    return labels

def Unnorm_HSC_LaplacianEigenmaps(A, n_clusters = None, plot_eig=False, plot_fiedler=False, norm=False, embedding="LaplacianEigenmaps", d=2, beta=None):
    labels = hardSpectralClustering(A, n_clusters, plot_eig, plot_fiedler, norm, embedding, d, beta)
    return labels
    
def Norm_HSC_LaplacianEigenmaps(A, n_clusters = None, plot_eig=False, plot_fiedler=False, norm=True, embedding="LaplacianEigenmaps", d=2, beta=None):
    labels = hardSpectralClustering(A, n_clusters, plot_eig, plot_fiedler, norm, embedding, d, beta)
    return labels
    
def Unnorm_SSC_LaplacianEigenmaps(A, n_clusters = None, plot_eig=False, plot_fiedler=False, norm=False, embedding="LaplacianEigenmaps", d=2, beta=None):
    labels = softSpectralClustering(A, n_clusters, plot_eig, plot_fiedler, norm, embedding, d, beta)
    return labels
    
def Norm_SSC_LaplacianEigenmaps(A, n_clusters = None, plot_eig=False, plot_fiedler=False, norm=True, embedding="LaplacianEigenmaps", d=2, beta=None):
    labels = softSpectralClustering(A, n_clusters, plot_eig, plot_fiedler, norm, embedding, d, beta)
    return labels
    
def HSC_ShiMalik_LaplacianEigenmaps(A, n_clusters = None, plot_eig=False, plot_fiedler=False, embedding="LaplacianEigenmaps", d=2, beta=None):
    labels = HSC_ShiMalik(A, n_clusters, plot_eig, plot_fiedler, embedding, d, beta)
    return labels
    
def SSC_ShiMalik_LaplacianEigenmaps(A, n_clusters = None, plot_eig=False, plot_fiedler=False, embedding="LaplacianEigenmaps", d=2, beta=None):
    labels = SSC_ShiMalik(A, n_clusters, plot_eig, plot_fiedler, embedding, d, beta)
    return labels


