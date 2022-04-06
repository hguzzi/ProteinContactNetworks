import os
import networkx as nx
import pcn_final
import numpy as np
from scipy.linalg import eigh
import pcn_pymol_scripts
from sys import platform
import configparser
import tkinter as tk
from tkinter import ttk

def userAnalysisChoice():
    
    global previous_scene 
    global isback
    global back_button
    global reset_button
    global analysis_fr
    global comp_adj_fr
     
    previous_scene = "select_proteins"
    comp_adj_fr.pack_forget()
    
    if isback or isreset:
       
        analysis_fr = tk.Frame(window)
 
    proteins_to_analyze_tk = tk.Label(analysis_fr, text = "List of proteins to analyze: {}".format(str(proteins_list)))
    proteins_to_analyze_tk.pack() 
    method_analysis_tk = tk.Label(analysis_fr, text="Choose the method to use for the analysis of the PCNs: ")
    method_analysis_tk.pack()
    
    sc_button = tk.Button(
                        analysis_fr,
                        text="Spectral Clustering",
                        bg = button_bg,
                        fg = button_fg,
                        command = spectralClustering 
                         )
    sc_button.pack()                
            
    ec_button = tk.Button(
                        analysis_fr,
                        text="Embedding + Clustering",
                        bg = button_bg,
                        fg = button_fg,
                        command = embeddingClustering 
                         )
    ec_button.pack() 
         
    cd_button = tk.Button(
                        analysis_fr,
                        text="Community Detection",
                        bg = button_bg,
                        fg = button_fg,
                        command = communityDetection 
                         )
    cd_button.pack() 
         
    cm_button = tk.Button(
                        analysis_fr,
                        text="Centrality Analysis",
                        bg = button_bg,
                        fg = button_fg,
                        command = centralityMeasure
                         )
          
    cm_button.pack()

    back_button = tk.Button(analysis_fr, text='Back', command = back, bg = button_bg, fg = button_fg)
    reset_button = tk.Button(analysis_fr, text='Reset Analysis', command = reset, bg = button_bg, fg = button_fg) 
    back_button.pack(side=tk.RIGHT)
    reset_button.pack(side=tk.LEFT)         
      
    analysis_fr.pack()

def compute_adjs():
    
    global min_
    global max_
    global proteins_list
    global adj_matrixs
    global proteins_residue_names
    global proteins_tk
    global choice_fr
    global comp_adj_fr
    global protein_path
    global adj_filespath
    global window 
    
    proteins = str(proteins_tk.get())
    [proteins_list.append(protein.casefold()) for protein in proteins.replace(" ","").split(",") if protein.casefold() not in proteins_list] #strip space
    min_ = float(str(min_tk.get()))
    max_ = float(str(max_tk.get()))
    
    choice_fr.pack_forget()  
    window.update()
    
    pdb_list = proteins_list.copy()
    for i, protein in enumerate(pdb_list):
        if (not protein.endswith('.pdb')):
            pdb_list[i] = protein+'.pdb'
    
    pcn_final.checkIfFilesExists(pdb_list, "pdb", proteins_path)   
    
    for protein in proteins_list:
                
        computing_A_label = tk.Label(comp_adj_fr, text="computing adjacency matrix for protein {}... (This may take time)".format(protein))    
        computing_A_label.pack()
        comp_adj_fr.pack()
        window.update()
        
        protein_path = proteins_path+protein+".pdb"
        atoms = pcn_final.readPDBFile(protein_path) 
        residues = pcn_final.getResidueDistance(atoms) 
        dict_residue_names = pcn_final.associateResidueName(residues)
        proteins_residue_names[protein] = np.array(list (dict_residue_names.items()))              
        
        adj_matrixs[protein] = pcn_final.adjacent_matrix(output_path, residues, protein, min_, max_, comp_adj_fr, window)
                
    userAnalysisChoice()
        

def read_adjs():
    
    global min_
    global max_
    global proteins_list
    global adj_matrixs
    global proteins_tk
    global choice_fr
    global comp_adj_fr
    global proteins_residue_names
    global proteins_path
    global adj_filespath
    global window 
    
    proteins = str(proteins_tk.get())
    [proteins_list.append(protein.casefold()) for protein in proteins.replace(" ","").split(",") if protein.casefold() not in proteins_list] #strip space
    min_ = float(str(min_tk.get()))
    max_ = float(str(max_tk.get()))
    
    choice_fr.pack_forget()
    window.update()
    
    adj_list = [protein.casefold()+"_adj_mat_{}_{}.txt".format(min_, max_) for protein in proteins_list]
    pcn_final.checkIfFilesExists(adj_list, "adj", proteins_path, adj_filespath)
    
    for protein in proteins_list:
        
        protein_path = proteins_path+protein+".pdb"
        atoms = pcn_final.readPDBFile(protein_path) 
        residues = pcn_final.getResidueDistance(atoms) 
        dict_residue_names = pcn_final.associateResidueName(residues)
        proteins_residue_names[protein] = np.array(list (dict_residue_names.items()))              
        
        reading_A_label = tk.Label(comp_adj_fr, text="reading adjacency matrix for protein {}... ".format(protein))
        reading_A_label.pack()
        comp_adj_fr.pack()
        window.update()
        
        adj_matrixs[protein] = pcn_final.read_adj_mat(adj_filespath, protein, min_, max_)

    userAnalysisChoice()

    
def centralityMeasure():
    
    global previous_scene
    global isback 
    global back_button
    global reset_button
    global cm_fr
    global mb_cm
    global analysis_fr
    global checked
    global button_bg 
    global button_fg
    global supported_centralities_measures
    
    previous_scene = "userAnalysisChoice"
    if isback or isreset:
        cm_fr = tk.Frame(window)
        mb_cm = tk.Menubutton(cm_fr, text="Click here and Choose at least one centrality measure algorithm to use: ", bg = button_bg, fg = button_fg)
    
    analysis_fr.pack_forget()
        
    mb_cm.grid()
    mb_cm.menu = tk.Menu(mb_cm)
    mb_cm["menu"] = mb_cm.menu

    for i, centralities_measure in enumerate(supported_centralities_measures):
            
        checked_tk = tk.BooleanVar()
        mb_cm.menu.add_checkbutton(label = centralities_measure, onvalue = True, offvalue = False, variable = checked_tk)
        checked.append(checked_tk)            
            
    run_button = tk.Button(cm_fr, text="Run", bg = button_bg, fg = button_fg, command=centralityMeasureRun)
    back_button = tk.Button(cm_fr, text='Back', bg = button_bg, fg = button_fg, command = back)
    reset_button = tk.Button(cm_fr, text='Reset Analysis', bg = button_bg, fg = button_fg, command = reset)
        
    cm_fr.pack()
    mb_cm.pack()
        
    run_button.pack()
    back_button.pack()
    reset_button.pack()   

def spectralClustering():
    
    global previous_scene
    global isback 
    global back_button
    global reset_button
    global sc_fr
    global parameters_fr
    global analysis_fr
    global mb_sc
    global checked
    global ks_tk
    global button_bg 
    global button_fg
    global supported_algorithms_clustering
    
    previous_scene = "userAnalysisChoice"
    if isback or isreset:
        sc_fr = tk.Frame(window)
        parameters_fr = tk.Frame(window)
        ks_tk = tk.Entry(parameters_fr, text="Number of clusters for clustering", width=100)
        mb_sc = tk.Menubutton(sc_fr, text="Click here and Choose at least one spectral clustering algorithm to use: ", bg = button_bg, fg = button_fg)
    
    analysis_fr.pack_forget()
        
    mb_sc.grid()
    mb_sc.menu = tk.Menu(mb_sc)
    mb_sc["menu"] = mb_sc.menu
      
    for i, algorithm_clustering in enumerate(supported_algorithms_clustering):
        checked_tk = tk.BooleanVar()
        mb_sc.menu.add_checkbutton(label = algorithm_clustering, onvalue = True, offvalue = False, variable = checked_tk)
        checked.append(checked_tk)            
          
    ks_insert_label = tk.Label(parameters_fr, text="Enter number of clusters for spectral clustering: Enter an int, a list of ints splitted with ',': ")
    ks_insert_label.pack()
    ks_tk.pack()
    parameters_fr.pack()

    run_button = tk.Button(sc_fr, text="Run", bg = button_bg, fg = button_fg, command=spectralClusteringRun)
    back_button = tk.Button(sc_fr, text='Back', bg = button_bg, fg = button_fg, command = back)
    reset_button = tk.Button(sc_fr, text='Reset Analysis', bg = button_bg, fg = button_fg, command = reset) 

    sc_fr.pack()
    mb_sc.pack()
        
    run_button.pack()
    back_button.pack()
    reset_button.pack()    
    
def communityDetection():
    
    global previous_scene
    global isback 
    global back_button
    global reset_button
    global cd_fr
    global parameters_fr
    global analysis_fr
    global mb_cd
    global checked
    global ks_tk    
    global button_bg 
    global button_fg
    global supported_algorithms_communities
    
    previous_scene = "userAnalysisChoice"
    if isback or isreset:

        cd_fr = tk.Frame(window)
        mb_cd = tk.Menubutton(cd_fr, text="Click here and Choose at least one community detection algorithm to use: ", bg = button_bg, fg = button_fg)
        parameters_fr = tk.Frame(window)
        ks_tk = tk.Entry(parameters_fr, text="Number of clusters for clustering", width=100)

    analysis_fr.pack_forget()
     
    mb_cd.grid()
    mb_cd.menu = tk.Menu(mb_cd)
    mb_cd["menu"] = mb_cd.menu
      
    for i, algorithm_community in enumerate(supported_algorithms_communities):
        checked_tk = tk.BooleanVar()
        mb_cd.menu.add_checkbutton(label = algorithm_community, onvalue = True, offvalue = False, variable = checked_tk)
        checked.append(checked_tk)            
        
    ks_insert_label = tk.Label(parameters_fr, text="Enter number of communities for Async FluidC algorithm: Enter an int or a list of ints splitted with ',': ")
    ks_insert_label.pack()
    ks_tk.pack()
    parameters_fr.pack()
          
    run_button = tk.Button(cd_fr, text="Run", bg = button_bg, fg = button_fg, command=communityDetectionRun)
    back_button = tk.Button(cd_fr, text='Back', bg = button_bg, fg = button_fg, command = back)
    reset_button = tk.Button(cd_fr, text='Reset Analysis', bg = button_bg, fg = button_fg, command = reset)
     
    cd_fr.pack()
    mb_cd.pack()
     
    run_button.pack()
    back_button.pack()
    reset_button.pack()   

def embeddingClustering():

    global previous_scene
    global isback 
    global back_button
    global reset_button
    global ec_fr
    global parameters_fr
    global analysis_fr
    global mb_ec
    global ks_tk
    global d_tk
    global beta_tk
    global num_walks_tk
    global walk_len_tk
    global button_bg 
    global button_fg
    global supported_algorithms_embeddings
    
    previous_scene = "userAnalysisChoice"

    if isback or isreset:
    
        ec_fr = tk.Frame(window)
        parameters_fr = tk.Frame(window)
        ks_tk = tk.Entry(parameters_fr, text="Number of clusters for clustering", width=100)
        d_tk = tk.Entry(parameters_fr, text="d-dimensional embedding: ")
        beta_tk = tk.Entry(parameters_fr, text="decay factor HOPE embedding: ")
        num_walks_tk = tk.Entry(parameters_fr, text="Random walk lenght: ")
        walk_len_tk = tk.Entry(parameters_fr, text="Number of walks per node: ")
        mb_ec = tk.Menubutton(ec_fr, text="Click here and Choose at least one embedding+clustering algorithm to use: ", bg = button_bg, fg = button_fg)

    analysis_fr.pack_forget()
    
    mb_ec.grid()
    mb_ec.menu = tk.Menu(mb_ec)
    mb_ec["menu"] = mb_ec.menu

    for i, algorithm_emb_clustering in enumerate(supported_algorithms_embeddings):
        checked_tk = tk.BooleanVar()
        mb_ec.menu.add_checkbutton(label = algorithm_emb_clustering, onvalue = True, offvalue = False, variable = checked_tk)
        checked.append(checked_tk)  
                  
    ks_insert_label = tk.Label(parameters_fr, text="Enter number of clusters for embedding + clustering algorithm: Enter an int or a list of ints splitted with ',': ")
    ks_insert_label.pack()
    ks_tk.pack()
            
    d_insert_label = tk.Label(parameters_fr, text=" Enter d parameter for d-dimensional embedding: ")
    d_insert_label.pack()
    d_tk.pack()
       
    beta_insert_label = tk.Label(parameters_fr, text="Enter beta parameter for HOPE embedding: ")
    beta_insert_label.pack()
    beta_tk.pack()
            
    num_walks_insert_label = tk.Label(parameters_fr, text="Enter the lenght of each random walk: ")
    num_walks_insert_label.pack()
    num_walks_tk.pack()
     
    walk_len_insert_label = tk.Label(parameters_fr, text="Enter the number of walks per node: ")
    walk_len_insert_label.pack()
    walk_len_tk.pack()
    
    parameters_fr.pack()
        
    run_button = tk.Button(ec_fr, text="Run", bg = button_bg, fg = button_fg, command=embeddingClusteringRun)
    back_button = tk.Button(ec_fr, text='Back', bg = button_bg, fg = button_fg, command = back)
    reset_button = tk.Button(ec_fr, text='Reset Analysis', bg = button_bg, fg = button_fg, command = reset)
          
    ec_fr.pack()
    mb_ec.pack()
    run_button.pack()
    back_button.pack()
    reset_button.pack()   

def computeOrReadPCN():
    
    global choiceVar 
    
    if choiceVar.get() == 'adj':
        read_adjs()
    else: #=='pcn'
        compute_adjs()
        
def centralityMeasureRun():
    
    global previous_scene
    global checked
    global cm_fr 
    global adj_matrixs
    global proteins_residue_names
    global protein_path 
    global output_path
    global run_fr
    global window
    global supported_centralities_measures
    global add_slash_to_path
    
    cm_touse = []
    
    for i, check_button in enumerate(checked):
        if check_button.get(): #checked
            cm_touse.append(supported_centralities_measures[i])
            
    filepaths = []
    cm_fr.pack_forget()
    for protein, adj in adj_matrixs.items():
        
        protein_path = proteins_path+protein+".pdb"
        G = nx.from_numpy_matrix(adj)  
        residue_names = proteins_residue_names[protein]
        residue_names_1 = np.array(residue_names[:, 1], dtype = str)  
        run_fr.pack()
        for centrality_choice in cm_touse:
            
            method_to_call = getattr(pcn_final, centrality_choice) 

            compute_tk = tk.Label(run_fr, text="Compute {} node centrality for protein {}...".format(centrality_choice, protein)) 
            compute_tk.pack()
            window.update()
            centrality_measures = method_to_call(G, residue_names_1)#call the supported method from the pcn_final file
            
            pcn_final.save_centralities(output_path, centrality_measures, protein, centrality_choice) #save a txt file 

            plot_tk = tk.Label(run_fr, text="Plot {} node centrality for protein {}".format(centrality_choice, protein))
            plot_tk.pack()
            window.update()
            pcn_pymol_scripts.pymol_plot_centralities(output_path, centrality_measures, protein_path, centrality_choice, run_fr, window) #plot and save centralities with pymol
            filepath = "{}Centralities{}{}{}Sessions{}{}_{}_session.pse".format(output_path, add_slash_to_path, centrality_choice, add_slash_to_path, add_slash_to_path, protein, centrality_choice)
            filepaths.append(filepath)
            
    showResults(filepaths)
    
def spectralClusteringRun():
    
    global previous_scene
    global ks
    global ks_tk
    global checked
    global sc_fr 
    global adj_matrixs
    global proteins_residue_names
    global protein_path 
    global output_path
    global run_fr
    global window 
    global supported_algorithms_clustering
    global add_slash_to_path
    
    sc_touse = []
    
    for i, check_button in enumerate(checked):
        if check_button.get(): #checked
            sc_touse.append(supported_algorithms_clustering[i])        
    sc_fr.pack_forget()
    
    k_choice = str(ks_tk.get())
    if k_choice == 'best_k':
        #use max eigengap method to compute the best k to use
        pass
    elif k_choice.split(','):    
        ks = [int(k) for k in k_choice.replace(" ","").split(",")] #strip space
  
    filepaths = []
    parameters_fr.pack_forget()
    run_fr.pack()
    for protein, adj in adj_matrixs.items():
        G = nx.from_numpy_matrix(adj)
        protein_path = proteins_path+protein+".pdb"
        residue_names = proteins_residue_names[protein]
        residue_names_1 = np.array(residue_names[:, 1], dtype = str)
        for algorithm_choice in sc_touse:
          
            for k in ks:
                
                method_to_call = getattr(pcn_final, algorithm_choice)
                
                compute_tk = tk.Label(run_fr, text="Compute {} spectral clustering with k = {} for protein {}...".format(algorithm_choice, k, protein)) 
                compute_tk.pack()
                window.update()
                labels = method_to_call(adj, n_clusters=k)
                pcn_final.save_labels(output_path, labels, residue_names, protein, algorithm_choice, d, beta, walk_len, num_walks)
                
                plot_tk = tk.Label(run_fr, text="Plot {} spectral clustering with k = {} for protein {}...".format(algorithm_choice, k, protein))
                plot_tk.pack()
                window.update()
                pcn_pymol_scripts.pymol_plot(protein_path, output_path, "Clusters", algorithm_choice, k, run_fr, window)
                filepath = ""
                filepaths.append(filepath)
                
                plot_p_tk = tk.Label(run_fr, text="Compute and Plot partecipation coefficients with {} spectral clustering and k = {} for protein {}...".format(algorithm_choice, k, protein))
                plot_p_tk.pack()
                window.update()
                p = pcn_final.participation_coefs(G, labels, residue_names_1)
                pcn_final.save_part_coef(output_path, p, protein, algorithm_choice, k)
                output_path_p = "{}{}{}{}".format(output_path, add_slash_to_path, algorithm_choice, add_slash_to_path)
                pcn_pymol_scripts.pymol_plot_part_coefs(p, protein_path, output_path_p, algorithm_choice, k, run_fr, window)
                filepath_p = "{}Part_coefs_Sessions{}{}_part_coefs_{}_k{}_session.pse".format(output_path_p, add_slash_to_path, protein, algorithm_choice, k)
                filepaths.append(filepath_p)
                
    showResults(filepaths)
    
def communityDetectionRun():
    
    global previous_scene
    global ks 
    global ks_tk
    global checked
    global cd_fr 
    global adj_matrixs
    global proteins_residue_names
    global protein_path 
    global output_path
    global run_fr
    global window
    global supported_algorithms_communities
    global add_slash_to_path
    
    cd_touse = []
    for i, check_button in enumerate(checked):
        if check_button.get(): #checked
            cd_touse.append(supported_algorithms_communities[i])
    cd_fr.pack_forget()
    
    if "asyn_fluidc" in cd_touse:
        k_choice = str(ks_tk.get())
        if k_choice.split(','):    
            ks = [int(k) for k in k_choice.replace(" ","").split(",")] #strip space
    
    parameters_fr.pack_forget()
    
    filepaths = []
    for protein, adj in adj_matrixs.items():
        protein_path = proteins_path+protein+".pdb"
        G = nx.from_numpy_matrix(adj) 
        residue_names = proteins_residue_names[protein]
        residue_names_1 = np.array(residue_names[:, 1], dtype = str)
        
        run_fr.pack()
        for algorithm_choice in cd_touse:
            method_to_call = getattr(pcn_final, algorithm_choice)
                   
            #if the algorithm is asyn fluidic
            if (algorithm_choice == 'asyn_fluidc'):
                for k in ks:
                            
                    compute_tk = tk.Label(run_fr, text="Compute Asyn FluidC communities with k = {} for protein {}...".format(k, protein)) 
                    compute_tk.pack()
                    window.update()
                    labels = method_to_call(G, k) #call the method
                    pcn_final.save_labels(output_path, labels, residue_names, protein,  method=algorithm_choice) #save the communities as txt file
                    
                    plot_tk = tk.Label(run_fr, text="Plot Asyn FluidC with k = {} for protein {}...".format(k, protein))
                    plot_tk.pack()
                    window.update()
                    pcn_pymol_scripts.pymol_plot(protein_path, output_path, "Communities", algorithm_choice, k, run_fr, window) #plot and save the communities with pymol
                    filepath = "{}{}{}Sessions{}{}_{}_{}_{}{}_session.pse".format(output_path, algorithm_choice, add_slash_to_path, add_slash_to_path, protein, "Communities", algorithm_choice, "ncoms", k)
                    filepaths.append(filepath)
                    plot_p_tk = tk.Label(run_fr, text="Compute and Plot partecipation coefficients with Asyn FluidC and k = {} for protein {}...".format(k, protein))
                    plot_p_tk.pack()
                    window.update()                   
                    p = pcn_final.participation_coefs(G, labels, residue_names_1)
                    pcn_final.save_part_coef(output_path, p, protein, algorithm_choice, k) #save the part coefs as txt file
                    output_path_p = "{}{}{}".format(output_path, algorithm_choice, add_slash_to_path)
                    pcn_pymol_scripts.pymol_plot_part_coefs(p, protein_path, output_path_p, algorithm_choice, k, run_fr, window) #plot and save part coefs with pymol
                    filepath_p = "{}Part_coefs_Sessions{}{}_part_coefs_{}_k{}_session.pse".format(output_path_p, add_slash_to_path, protein, algorithm_choice, k)
                    filepaths.append(filepath_p)
            
            else:#if the community detection algorithm is not Asyn Fluidc, no need to specify the number of communities
                
                compute_tk = tk.Label(run_fr, text="Compute {} communities for protein {}...".format(algorithm_choice, protein)) 
                compute_tk.pack()
                window.update()
                labels = method_to_call(G) #call the method 
                n_coms = int( max(labels) + 1)
                pcn_final.save_labels(output_path, labels, residue_names, protein,  method=algorithm_choice) #save communities as txt 
                
                plot_tk = tk.Label(run_fr, text="Plot {} with ncoms = {} for protein {}...".format(algorithm_choice, n_coms, protein))
                plot_tk.pack()
                window.update()
                pcn_pymol_scripts.pymol_plot(protein_path, output_path, "Communities", algorithm_choice, n_coms, run_fr, window) #plot and save communities with pymol
                filepath = "{}{}{}Sessions{}{}_{}_{}_{}{}_session.pse".format(output_path, algorithm_choice, add_slash_to_path, add_slash_to_path, protein, "Communities", algorithm_choice, "ncoms", n_coms)
                filepaths.append(filepath)
                
                plot_p_tk = tk.Label(run_fr, text="Compute and Plot partecipation coefficients with {} and ncoms = {} for protein {}...".format(algorithm_choice, n_coms, protein))
                plot_p_tk.pack()
                window.update()                   
                p = pcn_final.participation_coefs(G, labels, residue_names_1)
                pcn_final.save_part_coef(output_path, p, protein, algorithm_choice, n_coms)
                output_path_p = "{}{}{}".format(output_path, algorithm_choice, add_slash_to_path)
                pcn_pymol_scripts.pymol_plot_part_coefs(p, protein_path, output_path_p, algorithm_choice, n_coms, run_fr, window)
                filepath_p = "{}Part_coefs_Sessions{}{}_part_coefs_{}_k{}_session.pse".format(output_path_p, add_slash_to_path, protein, algorithm_choice, n_coms)
                filepaths.append(filepath_p)
      
    showResults(filepaths)
    
def embeddingClusteringRun():
    
    global previous_scene
    global ks 
    global ks_tk
    global d_tk
    global d
    global beta
    global beta_tk
    global num_walks
    global num_walks_tk
    global walk_len
    global walk_len_tk
    global checked
    global ec_fr 
    global adj_matrixs
    global run_fr 
    global window
    global proteins_residue_names
    global protein_path 
    global output_path
    global supported_algorithms_embeddings
    global add_slash_to_path
    
    ec_touse = []
    for i, check_button in enumerate(checked):
        if check_button.get(): #checked
            ec_touse.append(supported_algorithms_embeddings[i])
    ec_fr.pack_forget()
    
    k_choice = str(ks_tk.get())
    if k_choice.split(','):    
        ks = [int(k) for k in k_choice.replace(" ","").split(",")] #strip space
    
    d = int(d_tk.get())
            
    filepaths = []
    parameters_fr.pack_forget()
    for protein, adj in adj_matrixs.items():
        protein_path = proteins_path+protein+".pdb"
        G = nx.from_numpy_matrix(adj) 
        residue_names = proteins_residue_names[protein]
        residue_names_1 = np.array(residue_names[:, 1], dtype = str)
        for algorithm_choice in ec_touse:
           
            if "hope" in algorithm_choice:
                beta = float(beta_tk.get())
            if "node2vec" in algorithm_choice:
                num_walks = int(num_walks_tk.get())
                walk_len =  int(walk_len_tk.get())
                
            for k in ks:

                method_to_call = getattr(pcn_final, algorithm_choice)
                compute_tk = tk.Label(run_fr, text="Compute {} embedding + clustering with k = {}, d = {},  beta = {}, num_walks = {} and walk_len = {} for protein {}...".format(algorithm_choice, k, d, beta, num_walks, walk_len, protein)) 
                compute_tk.pack()
                run_fr.pack()
                window.update()
                labels = method_to_call(adj, n_clusters=k, d=d, beta=beta, walk_len=walk_len, num_walks=num_walks)
                pcn_final.save_labels(output_path, labels, residue_names, protein, algorithm_choice, d, beta, walk_len, num_walks)
                
                plot_tk = tk.Label(run_fr, text="Plot {} embedding + clustering with k = {}, d = {},  beta = {}, num_walks = {} and walk_len = {} for protein {}...".format(algorithm_choice, k, d, beta, num_walks, walk_len, protein)) 
                plot_tk.pack()
                run_fr.pack()
                window.update()
                pcn_pymol_scripts.pymol_plot_embeddings(protein_path, output_path, "ClustersEmbeddings", algorithm_choice, k, d, beta, walk_len, num_walks, run_fr, window)
                                
                if "hope" in algorithm_choice:
                    filepath = "{}{}{}Sessions{}{}_{}_d{}_beta{}_k{}_session.pse".format(output_path, algorithm_choice, add_slash_to_path, add_slash_to_path, protein, algorithm_choice, d, beta, k)
                    
                elif "node2vec" in algorithm_choice:
                    filepath = "{}{}{}Sessions{}{}_{}_d{}_wl{}_nw{}_k{}_session.pse".format(output_path, algorithm_choice, add_slash_to_path, add_slash_to_path, protein, algorithm_choice, d, walk_len, num_walks, k)

                else:
                    filepath = "{}{}{}Sessions{}{}_{}_d{}_k{}_session.pse".format(output_path, algorithm_choice, add_slash_to_path, add_slash_to_path, protein, algorithm_choice, d, k)
                    
                filepaths.append(filepath)
                
                plot_p_tk = tk.Label(run_fr, text="Compute and Plot partecipation coefficients with {}, k = {}, d = {},  beta = {}, num_walks = {} and walk_len = {} for protein {}...".format(algorithm_choice, k, d, beta, num_walks, walk_len, protein)) 
                plot_p_tk.pack()
                run_fr.pack()
                window.update() 
                p = pcn_final.participation_coefs(G, labels, residue_names_1)
                pcn_final.save_part_coef(output_path, p, protein, algorithm_choice, k)
                output_path_p = "{}{}{}{}".format(output_path, add_slash_to_path, algorithm_choice, add_slash_to_path)
                pcn_pymol_scripts.pymol_plot_part_coefs(p, protein_path, output_path_p, algorithm_choice, k, run_fr, window)
                filepath_p = "{}Part_coefs_Sessions{}{}_part_coefs_{}_k{}_session.pse".format(output_path_p, add_slash_to_path, protein, algorithm_choice, k)
                filepaths.append(filepath_p)
    
    showResults(filepaths)

def reset():
   
    global isreset
    global checked
    global ks
    global d
    global beta
    global num_walks
    global walk_len
        
    isreset = True
    
    comp_adj_fr.pack_forget()
    analysis_fr.pack_forget()
    cm_fr.pack_forget()
    cd_fr.pack_forget()
    sc_fr.pack_forget()
    ec_fr.pack_forget()
    parameters_fr.pack_forget()
    run_fr.pack_forget()
    results_fr.pack_forget()
    
    checked = list()
    ks = list()
    d = None
    beta = None
    num_walks = None
    walk_len = None
    
    select_proteins()
    
def back():
    
    global isback
    global previous_scene
    global comp_adj_fr
    global analysis_fr
    global cm_fr
    global cd_fr
    global sc_fr
    global ec_fr
    global parameters_fr
    
    isback = True 
    if previous_scene == "select_proteins":
        comp_adj_fr.pack_forget()
        analysis_fr.pack_forget()
        select_proteins()
    if previous_scene == "userAnalysisChoice":
        cm_fr.pack_forget()
        cd_fr.pack_forget()
        sc_fr.pack_forget()
        ec_fr.pack_forget()
        parameters_fr.pack_forget()
        userAnalysisChoice()

def showResults(filepaths):
    
    global run_fr
    global results_fr
    global window 
    global reset_button
    run_fr.pack_forget()
    
    for filepath in filepaths:
       
       filename = os.path.basename(filepath)
       
       label = tk.Label(results_fr, text = filename)
       label.pack()
       button = tk.Button(results_fr, text = "Open PyMOL session", bg = button_bg, fg = button_fg, command = os.system(filepath))
       button.pack()
       reset_button = tk.Button(results_fr, text = "Start a new analysis", bg = button_bg, fg = button_fg, command = reset)
       results_fr.pack()    
       window.update()
    
    reset_button.pack(side=tk.BOTTOM)
    results_fr.pack()    
    window.update()
    
    
def select_proteins():
     
    global isreset
    global isback 
    global output_path
    global proteins_tk
    global min_tk
    global max_tk
    global proteins_path
    global adj_filespath
    global choice_fr
    global initial_fr_button
    global comp_adj_fr
    
    if isback:
        
        choice_fr = tk.Frame(window)
        comp_adj_fr = tk.Frame(window)
        proteins_tk = tk.Entry(choice_fr, text='PDBs codes:', width=100)
        min_tk = tk.Entry(choice_fr, text='Min threshold:')
        max_tk = tk.Entry(choice_fr, text='Max threshold:')
        
    if isreset:    
    
        global analysis_fr
        global cm_fr
        global sc_fr
        global cd_fr
        global ec_fr
        global parameters_fr
        global run_fr
        global results_fr
        
        choice_fr = tk.Frame(window)
        analysis_fr = tk.Frame(window)
        cm_fr = tk.Frame(window)
        sc_fr = tk.Frame(window)
        cd_fr = tk.Frame(window)
        ec_fr = tk.Frame(window)
        parameters_fr = tk.Frame(window)
        comp_adj_fr = tk.Frame(window)
        run_fr = tk.Frame(window)
        results_fr = tk.Frame(window)
    
        proteins_tk = tk.Entry(choice_fr, text='PDBs codes:', width=100)
        min_tk = tk.Entry(choice_fr, text='Min threshold:')
        max_tk = tk.Entry(choice_fr, text='Max threshold:')

    initial_fr_button.pack_forget()
    if(os.path.isfile(os.getcwd()+add_slash_to_path+"config.ini")):
          
        config_obj = configparser.ConfigParser()
        config_obj.read(os.getcwd()+add_slash_to_path+"config.ini")
        paths = config_obj["user_paths"]
        output_path = paths["output_path"]
        proteins_path = paths["proteins_path"]
        adj_filespath = paths["adj_filespath"]
        paths = tk.Label(choice_fr, text="Paths in the config file: \n output path = {} \n proteins_path = {} \n adj_filespath = {}".format(output_path, proteins_path, adj_filespath))
        paths.pack(pady=0)
      
    else:
           
        proteins_path_tk = tk.Entry(choice_fr, text='PDB Directory:')
        proteins_path_tk.pack()
        proteins_path = str(proteins_path_tk.get())
        adj_filespath_tk = tk.Entry(choice_fr, text='Adj Directory:')
        adj_filespath_tk.pack()
        adj_filespath = str(adj_filespath_tk.get())
        output_path_tk = tk.Entry(choice_fr, text='PDB Directory:')
        output_path_tk.pack()
        output_path = str(output_path_tk.get())
        
    proteins_path = proteins_path+add_slash_to_path
    output_path = output_path+add_slash_to_path
    adj_filespath = adj_filespath+add_slash_to_path
    proteins_insert_label = tk.Label(choice_fr, text="Please Insert Protein PDB Identifiers, separated by comma, without .pdb, e.g. 7nxc for 7nxc.pdb")
    proteins_insert_label.pack()
    proteins_tk.pack()
      
    adj_pdb_choice = tk.Label(choice_fr, text="Format File Input: PDB structures or Preprocessed PCN")
    adj_pdb_choice.pack()
     
    choices = ("pdb", "adj")
    choiceVar.set(choices[0])
    cb = ttk.Combobox(choice_fr, textvariable=choiceVar, values=choices)
    cb.pack()
      
    min_insert_label = tk.Label(choice_fr, text="Please enter non covalent bonds threshold distance for PCN costruction")
    min_insert_label.pack()  
    min_tk.pack()

    max_insert_label = tk.Label(choice_fr, text="Please enter only significant bonds threshold distance for PCN costruction")
    max_insert_label.pack()
    max_tk.pack()
      
    submit_button = tk.Button(
                            choice_fr,
                            text="Start compute PCNs",
                            bg = button_bg, 
                            fg = button_fg,
                            command = computeOrReadPCN
                             )
    submit_button.pack(pady=12)
    choice_fr.pack()

def main():
    
    global supported_algorithms_clustering
    supported_algorithms_clustering = ["unnorm_ssc", "norm_ssc", "unnorm_hsc", "norm_hsc", "hsc_shimalik", "ssc_shimalik", "skl_spectral_clustering"]
    global supported_algorithms_embeddings
    supported_algorithms_embeddings = [
                                       "unnorm_ssc_hope", "norm_ssc_hope", "unnorm_hsc_hope", "norm_hsc_hope", "hsc_shimalik_hope", "ssc_shimalik_hope",
                                       
                                       "unnorm_ssc_laplacianeigenmaps", "norm_ssc_laplacianeigenmaps", "unnorm_hsc_laplacianeigenmaps", "norm_hsc_laplacianeigenmaps",
                                       "hsc_shimalik_laplacianeigenmaps", "ssc_shimalik_laplacianeigenmaps", 
                                       
                                       "unnorm_ssc_node2vec", "norm_ssc_node2vec", "unnorm_hsc_node2vec", "norm_hsc_node2vec", "hsc_shimalik_node2vec", 
                                       "ssc_shimalik_node2vec"
                                       ]
    global supported_algorithms_communities
    supported_algorithms_communities = ["louvain", "leiden", "walktrap", "asyn_fluidc", "greedy_modularity", "infomap", "spinglass"]
    global supported_centralities_measures
    supported_centralities_measures = ["closeness", "eigenvector_c", "betweenness", "degree_c"]
    
    global add_slash_to_path
    if platform == "linux" or platform == "linux2":
        # linux
        add_slash_to_path = '/'                    #add_slash_to_path is path separator 
        os.system('clear')
    elif platform == "darwin":
        # OS X
        add_slash_to_path = '/'
        os.system('clear')
    elif platform == "win32":
        # Windows...
        add_slash_to_path = '\\'
        supported_algorithms_communities.remove("infomap")
    
    global output_path
    output_path = ""
    global proteins_path
    proteins_path = ""
    global adj_filespath
    adj_filespath = ""
    global proteins_list
    proteins_list = []
    global proteins
    proteins = ""
    global adj_matrixs
    adj_matrixs = dict() # {pdb: adj} 
    global proteins_residue_names
    proteins_residue_names = dict() #{pdb: residue_names}
    global min_
    min_ = 4.0
    global max_
    max_ = 8.0
    global previous_scene
    previous_scene = "select_proteins"
    global isback
    isback = False
    global isreset
    isreset = False
    
    global window
    window = tk.Tk()
    window.title("PCN-Miner 0.9.1-alpha")
    window.rowconfigure(0, minsize=800, weight=1)
    window.columnconfigure(1, minsize=800, weight=1)
    window.geometry('600x600+50+50')

    label = tk.Label(window, text = 'PCN-Miner 0.9.1-alpha', fg = 'black',
                     font = (None, 30), height = 2)
    label.pack(side = tk.TOP)
    #window.iconbitmap('./assets/pythontutorial.ico') #icon

    welcome = tk.Label(window, text="Protein Contact Network Miner 0.9.1-alpha \n Software Available under CC-BY Licence \n Free for Academic Usage ",
                       foreground="black")  
    welcome.pack()
    
    global choiceVar
    choiceVar = tk.StringVar()
    
    global choice_fr
    choice_fr = tk.Frame(window)
    global analysis_fr
    analysis_fr = tk.Frame(window)
    global cm_fr
    cm_fr = tk.Frame(window)
    global sc_fr
    sc_fr = tk.Frame(window)
    global cd_fr
    cd_fr = tk.Frame(window)
    global ec_fr
    ec_fr = tk.Frame(window)
    global comp_adj_fr
    comp_adj_fr = tk.Frame(window)
    global parameters_fr
    parameters_fr = tk.Frame(window)
    global run_fr
    run_fr = tk.Frame(window)
    global results_fr
    results_fr = tk.Frame(window)

    global proteins_tk
    proteins_tk = tk.Entry(choice_fr, text='PDBs codes:', width=100)
    global min_tk
    min_tk = tk.Entry(choice_fr, text='Min threshold:')
    global max_tk
    max_tk = tk.Entry(choice_fr, text='Max threshold:')
    
    global ks_tk
    ks_tk = tk.Entry(parameters_fr, text="Number of clusters for clustering", width=100)
    global d_tk
    d_tk = tk.Entry(parameters_fr, text="d-dimensional embedding: ")
    global beta_tk
    beta_tk = tk.Entry(parameters_fr, text="decay factor HOPE embedding: ")
    global num_walks_tk
    num_walks_tk = tk.Entry(parameters_fr, text="Random walk lenght: ")
    global walk_len_tk
    walk_len_tk = tk.Entry(parameters_fr, text="Number of walks per node: ")
    
    global button_bg
    button_bg = "cyan"
    global button_fg
    button_fg = "black"

    global mb_cm
    mb_cm = tk.Menubutton(cm_fr, text="Click here and Choose at least one centrality measure algorithm to use: ", bg = button_bg, fg = button_fg)
    global mb_ec
    mb_ec = tk.Menubutton(ec_fr, text="Click here and Choose at least one embedding+clustering algorithm to use: ", bg = button_bg, fg = button_fg)
    global mb_sc
    mb_sc = tk.Menubutton(sc_fr, text="Click here and Choose at least one spectral clustering algorithm to use: ", bg = button_bg, fg = button_fg)
    global mb_cd
    mb_cd = tk.Menubutton(cd_fr, text="Click here and Choose at least one community detection algorithm to use: ", bg = button_bg, fg = button_fg)

    global checked
    checked = list()
    
    global ks
    ks = list()
    global d
    d = None
    global beta
    beta = None
    global num_walks
    num_walks = None
    global walk_len
    walk_len = None
    
    global back_button
    back_button = None
    global reset_button
    reset_button = None

    global initial_fr_button
    initial_fr_button = tk.Frame(window)
    
    initial_button = tk.Button(
                                initial_fr_button,
                                text="Press this button to continue",
                                bg = button_bg,
                                fg = button_fg,
                                command = select_proteins
                              )

    initial_button.pack()
    initial_fr_button.pack()

    window.mainloop()

if __name__ == '__main__':
    main()