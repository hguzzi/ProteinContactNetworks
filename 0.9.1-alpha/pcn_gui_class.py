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


class PCNMinerGUI():
    """
    Graphic User Interface Object for PCN-Miner.
    """
    def __init__(self, master):
        """
        Initialize all attributes of the object
        """
        
        #initialize supported algorithms
        self.supported_algorithms_clustering = ["unnorm_ssc", "norm_ssc", "unnorm_hsc", "norm_hsc", "hsc_shimalik", "ssc_shimalik", "skl_spectral_clustering"]
        self.supported_algorithms_embeddings = [
                                                "fuzzycmeans_hope", "kmeans_hope", "fuzzycmeans_laplacianeigenmaps", "kmeans_laplacianeigenmaps" ,
                                                "fuzzycmeans_node2vec", "kmeans_node2vec"
                                                ]
        self.supported_algorithms_communities = ["louvain", "leiden", "walktrap", "asyn_fluidc", "greedy_modularity", "infomap", "spinglass"]
        self.supported_centralities_measures = ["closeness", "eigenvector_c", "betweenness", "degree_c"]
    
        #handle different path separator between different Operative Systems
        if platform == "linux" or platform == "linux2":
            # linux
            self.add_slash_to_path = '/'                    #add_slash_to_path is path separator 
            os.system('clear')
        elif platform == "darwin":
            # OS X
            self.add_slash_to_path = '/'
            os.system('clear')
        elif platform == "win32":
            # Windows...
            self.add_slash_to_path = '\\'
            self.supported_algorithms_communities.remove("infomap")
        
        #initialize paths
        self.output_path = ""
        self.proteins_path = ""
        self.adj_filespath = ""
        
        #initialize proteins
        self.proteins_list = []
        self.proteins = ""
        
        #dictionary {pdb: adj}, each protein (pdb) is associated with the adj matrix that represents the Protein Contact Network of the protein. 
        self.adj_matrixs = dict() 
        #dictionary {pdb: residue_names}, each protein (pdb) is associated with the list of the name of the residues of the protein 
        self.proteins_residue_names = dict() 
        
        #initialize min and max distance threshold for non covalent and significative interactions between amino acids 
        self.min_ = 4.0
        self.max_ = 8.0
        
        #save the previous "scene" seen by the user, is used when the user click the Back Button
        self.previous_scene = "select_proteins"
        
        #true when user click the Back button
        self.isback = False
        #true when user click the Reset button
        self.isreset = False
        #true if config file with directory paths exists
        self.config = True
                    
        #Initialize and set some parameters for the GUI window
        self.window = master
        self.window.title("PCN-Miner 0.9.1-alpha")
        self.window.rowconfigure(0, minsize=800, weight=1)
        self.window.columnconfigure(1, minsize=800, weight=1)
        self.window.geometry('600x600+50+50')

        #name of the program to display always on the top of the window
        self.label_program_name = tk.Label(self.window, text = 'PCN-Miner 0.9.1-alpha', fg = 'black',
                         font = (None, 30), height = 2)
        self.label_program_name.pack(side = tk.TOP)
        #window.iconbitmap('./assets/pythontutorial.ico') #icon

        #welcome message
        welcome = tk.Label(self.window, text="Protein Contact Network Miner 0.9.1-alpha \n Software Available under CC-BY Licence \n Free for Academic Usage ",
                           foreground="black")  
        welcome.pack()
        
        #string variable, "pdb" or "adj"
        self.choiceVar = tk.StringVar()
        
        #initialize frames
        self.choice_fr = tk.Frame(self.window)         #for choose proteins to analyze, min and max distance threshold
        self.analysis_fr = tk.Frame(self.window)       #for choose analysis approach on the PCN
        self.cm_fr = tk.Frame(self.window)             #for choose centralities measure algorithm
        self.sc_fr = tk.Frame(self.window)             #for choose spectral clustering algorithm   
        self.cd_fr = tk.Frame(self.window)             #for choose community detection algorithm
        self.ec_fr = tk.Frame(self.window)             #for choose embedding + clustering algorithm
        self.comp_adj_fr = tk.Frame(self.window)       #for compute the adjacency matrixs (PCNs)
        self.parameters_fr = tk.Frame(self.window)     #frame for entry parameters
        self.run_fr = tk.Frame(self.window)            #frame for run analysis on the PCN
        self.results_fr = tk.Frame(self.window)        #for display results
        
        #initialize entries
        
        #paths user entry, to use when no config file exist
        self.proteins_path_tk = tk.Entry(self.choice_fr, text='PDB Directory:', width=80)
        self.adj_filespath_tk = tk.Entry(self.choice_fr, text='Adj Directory:', width=80)
        self.output_path_tk = tk.Entry(self.choice_fr, text="Output Directory:", width=80)
        
        #entry proteins to analyze
        self.proteins_tk = tk.Entry(self.choice_fr, text='PDBs codes:', width=100)
        #entry min distance threshold
        self.min_tk = tk.Entry(self.choice_fr, text='Min threshold:')
        #entry max distance threshold
        self.max_tk = tk.Entry(self.choice_fr, text='Max threshold:')  
        #entry list of number of clusters to use
        self.ks_tk = tk.Entry(self.parameters_fr, text="Number of clusters for clustering", width=100)
        #entry parameter d for d-dimensional embedding
        self.d_tk = tk.Entry(self.parameters_fr, text="d-dimensional embedding: ")
        #entry parameter beta if the user wants to use HOPE embedding
        self.beta_tk = tk.Entry(self.parameters_fr, text="decay factor HOPE embedding: ")
        #entry parameters number of walks and walk lenght if the user wants to use node2vec embedding 
        self.num_walks_tk = tk.Entry(self.parameters_fr, text="Random walk lenght: ")
        self.walk_len_tk = tk.Entry(self.parameters_fr, text="Number of walks per node: ")
        
        #background and foreground button colors
        self.button_bg = "cyan"
        self.button_fg = "black"
        
        #initialize menu buttons that enables the user to select multiple algorithms to use in the analysis
        self.mb_cm = tk.Menubutton(self.cm_fr, text="Click here and Choose at least one centrality measure algorithm to use: ", bg = self.button_bg, fg = self.button_fg)
        self.mb_ec = tk.Menubutton(self.ec_fr, text="Click here and Choose at least one embedding+clustering algorithm to use: ", bg = self.button_bg, fg = self.button_fg)
        self.mb_sc = tk.Menubutton(self.sc_fr, text="Click here and Choose at least one spectral clustering algorithm to use: ", bg = self.button_bg, fg = self.button_fg)
        self.mb_cd = tk.Menubutton(self.cd_fr, text="Click here and Choose at least one community detection algorithm to use: ", bg = self.button_bg, fg = self.button_fg)
        
        #list of all selected algorithms  
        self.checked = list()
        
        #list of number of clusters given by the user
        self.ks = list()
        #parameter for d-dimensional embedding given by the user
        self.d = None
        #beta decay factor for HOPE embedding given by the user
        self.beta = None
        #number of walks per node for node2vec embedding given by the user
        self.num_walks = None
        #lenght of random walks for node2vec embedding given by the user
        self.walk_len = None
        
        #initialize back and reset button, the first "scene" ofcourse doesn't need a back or reset button
        self.back_button = None
        self.reset_button = None
        
        #Start Analysis Button
        self.initial_fr_button = tk.Frame(self.window)
            
        self.initial_button = tk.Button(
                                    self.initial_fr_button,
                                    text="Press this button to continue",
                                    bg = self.button_bg,
                                    fg = self.button_fg,
                                    command = self.select_proteins
                                  )

        self.initial_button.pack()
        self.initial_fr_button.pack()
            
    def userAnalysisChoice(self):
        """
            User have to choose one method between different alternatives to analyze the PCN: 
                Spectral Clustering;
                Community Detection;
                Embedding+Clustering;
                Centrality Measures. 
            
        """
        #if the user press back, return to the "select proteins" scene 
        self.previous_scene = "select_proteins"
        self.comp_adj_fr.pack_forget()
        
        if self.isback or self.isreset:
            #if this page is loaded because the user press the back or the reset button we have to re-initialize all the frames used in the function
            self.analysis_fr = tk.Frame(self.window)
        
        #Some text that helps the user 
        proteins_to_analyze_tk = tk.Label(self.analysis_fr, text = "List of proteins to analyze: {}".format(str(self.proteins_list)))
        proteins_to_analyze_tk.pack() 
        method_analysis_tk = tk.Label(self.analysis_fr, text="Choose the method to use for the analysis of the PCNs: ")
        method_analysis_tk.pack()
        
        #buttons that enables the user to select the method approach to use
        sc_button = tk.Button(
                            self.analysis_fr,
                            text="Spectral Clustering",
                            bg = self.button_bg,
                            fg = self.button_fg,
                            command = self.spectralClustering 
                             )
        sc_button.pack()                
                
        ec_button = tk.Button(
                            self.analysis_fr,
                            text="Embedding + Clustering",
                            bg = self.button_bg,
                            fg = self.button_fg,
                            command = self.embeddingClustering 
                             )
        ec_button.pack() 
             
        cd_button = tk.Button(
                            self.analysis_fr,
                            text="Community Detection",
                            bg = self.button_bg,
                            fg = self.button_fg,
                            command = self.communityDetection 
                             )
        cd_button.pack() 
             
        cm_button = tk.Button(
                            self.analysis_fr,
                            text="Centrality Analysis",
                            bg = self.button_bg,
                            fg = self.button_fg,
                            command = self.centralityMeasure
                             )
              
        cm_button.pack()
        
        #reset and back buttons
        self.back_button = tk.Button(self.analysis_fr, text='Back', command = self.back, bg = self.button_bg, fg = self.button_fg)
        self.reset_button = tk.Button(self.analysis_fr, text='Reset Analysis', command = self.reset, bg = self.button_bg, fg = self.button_fg) 
        self.back_button.pack(side=tk.RIGHT)
        self.reset_button.pack(side=tk.LEFT)         
          
        self.analysis_fr.pack()

    def computeOrReadPCN(self):
        
        if not self.config:
            self.proteins_path = str(self.proteins_path_tk.get())
            self.adj_filespath = str(self.adj_filespath_tk.get())
            self.output_path = str(self.output_path_tk.get())
            
            if (not self.proteins_path.endswith(self.add_slash_to_path)):
                self.proteins_path = self.proteins_path+self.add_slash_to_path
            if (not self.output_path.endswith(self.add_slash_to_path)):
                self.output_path = self.output_path+self.add_slash_to_path
            if (not self.adj_filespath.endswith(self.add_slash_to_path)):
                self.adj_filespath = self.adj_filespath+self.add_slash_to_path

        #if type file is adj -> read the adjs files
        if self.choiceVar.get() == 'adj':
            self.read_adjs()
        #if type file is pdb -> compute the adj matrix starting from the PDB file
        else: #=='pcn'
            self.compute_adjs()
    
    def compute_adjs(self):
        """
        Computes the adj matrix for all the proteins selected.
        """
        #populate the protein list from the user entry input
        proteins = str(self.proteins_tk.get())
        [self.proteins_list.append(protein.casefold()) for protein in proteins.replace(" ","").split(",") if protein.casefold() not in self.proteins_list] #strip space
        #save min and max distance threshold parameters 
        self.min_ = float(str(self.min_tk.get()))
        self.max_ = float(str(self.max_tk.get()))
        
        self.choice_fr.pack_forget()  
        self.window.update()
        
        #pdb codes list (adds .pdb to all the proteins)
        pdb_list = self.proteins_list.copy()
        for i, protein in enumerate(pdb_list):
            if (not protein.endswith('.pdb')):
                pdb_list[i] = protein+'.pdb'
        
        #check if the pdb files exists in the proteins directory, if a protein pdb doesn't exists in the directory the software download the pdb
        pcn_final.checkIfFilesExists(pdb_list, "pdb", self.proteins_path)   
        
        #for each protein in the protein list -> compute the adj matrix
        for protein in self.proteins_list:
                    
            computing_A_label = tk.Label(self.comp_adj_fr, text="computing adjacency matrix for protein {}... (This may take time)".format(protein))    
            computing_A_label.pack()
            self.comp_adj_fr.pack()
            self.window.update()
            
            protein_path = self.proteins_path+protein+".pdb"
            atoms = pcn_final.readPDBFile(protein_path) 
            residues = pcn_final.getResidueDistance(atoms) 
            dict_residue_names = pcn_final.associateResidueName(residues)
            self.proteins_residue_names[protein] = np.array(list (dict_residue_names.items()))              
            
            self.adj_matrixs[protein] = pcn_final.adjacent_matrix(self.output_path, residues, protein, self.min_, self.max_, self.comp_adj_fr, self.window)
        
        #next scene
        self.userAnalysisChoice()
            

    def read_adjs(self):
        """
        Read the adj matrix for all the proteins selected.
        """
        #populate the protein list from the user entry input
        proteins = str(self.proteins_tk.get())
        [self.proteins_list.append(protein.casefold()) for protein in proteins.replace(" ","").split(",") if protein.casefold() not in self.proteins_list] #strip space
        #save min and max distance threshold parameters 
        self.min_ = float(str(self.min_tk.get()))
        self.max_ = float(str(self.max_tk.get()))
        
        self.choice_fr.pack_forget()
        self.window.update()
        
        #list of adj matrix file name 
        adj_list = [protein.casefold()+"_adj_mat_{}_{}.txt".format(self.min_, self.max_) for protein in self.proteins_list]
        #check if the adj files exists in the proteins directory, if the adj of the protein doesn't exists in the directory the software will compute it
        pcn_final.checkIfFilesExists(adj_list, "adj", self.proteins_path, self.adj_filespath)
        
        #for each protein in the protein list -> read the adjacency matrix
        for protein in self.proteins_list:
            
            protein_path = self.proteins_path+protein+".pdb"
            atoms = pcn_final.readPDBFile(protein_path) 
            residues = pcn_final.getResidueDistance(atoms) 
            dict_residue_names = pcn_final.associateResidueName(residues)
            self.proteins_residue_names[protein] = np.array(list (dict_residue_names.items()))              
            
            reading_A_label = tk.Label(self.comp_adj_fr, text="reading adjacency matrix for protein {}... ".format(protein))
            reading_A_label.pack()
            self.comp_adj_fr.pack()
            self.window.update()
            
            self.adj_matrixs[protein] = pcn_final.read_adj_mat(self.adj_filespath, protein, self.min_, self.max_)
        
        #next scene
        self.userAnalysisChoice()

        
    def centralityMeasure(self):
        """
        Select the node centrality algorithms to use.
        """
        
        #if the user press back, return to the "userAnalysisChoice" scene 
        self.previous_scene = "userAnalysisChoice"
        
        if self.isback or self.isreset:
            #if this page is loaded because the user press the back or the reset button we have to re-initialize some attributes
            self.checked = list()
            self.cm_fr = tk.Frame(self.window)
            self.mb_cm = tk.Menubutton(self.cm_fr, text="Click here and Choose at least one centrality measure algorithm to use: ", bg = self.button_bg, fg = self.button_fg)
        
        self.analysis_fr.pack_forget()
        
        #configure the menu button
        self.mb_cm.grid()
        self.mb_cm.menu = tk.Menu(self.mb_cm)
        self.mb_cm["menu"] = self.mb_cm.menu
    
        #for each centrality algorithm in the list of supported centrality algorithms -> add a button in the menu
        for i, centralities_measure in enumerate(self.supported_centralities_measures):
                
            checked_tk = tk.BooleanVar()
            self.mb_cm.menu.add_checkbutton(label = centralities_measure, onvalue = True, offvalue = False, variable = checked_tk)
            self.checked.append(checked_tk)            
        
        #run, back and reset buttons
        run_button = tk.Button(self.cm_fr, text="Run", bg = self.button_bg, fg = self.button_fg, command=self.centralityMeasureRun)
        self.back_button = tk.Button(self.cm_fr, text='Back', bg = self.button_bg, fg = self.button_fg, command = self.back)
        self.reset_button = tk.Button(self.cm_fr, text='Reset Analysis', bg = self.button_bg, fg = self.button_fg, command = self.reset)
            
        self.cm_fr.pack()
        self.mb_cm.pack()
            
        run_button.pack()
        self.back_button.pack()
        self.reset_button.pack()   

    def spectralClustering(self):
        """
        Select the spectral clustering algorithms and numbers of clusters ks to use.
        """
        #if the user press back, return to the "userAnalysisChoice" scene 
        self.previous_scene = "userAnalysisChoice"
        if self.isback or self.isreset:
            #if this page is loaded because the user press the back or the reset button we have to re-initialize some attributes
            self.checked = list()
            self.sc_fr = tk.Frame(self.window)
            self.parameters_fr = tk.Frame(self.window)
            self.ks_tk = tk.Entry(self.parameters_fr, text="Number of clusters for clustering", width=100)
            self.mb_sc = tk.Menubutton(self.sc_fr, text="Click here and Choose at least one spectral clustering algorithm to use: ", bg = self.button_bg, fg = self.button_fg)
        
        self.analysis_fr.pack_forget()
            
        self.mb_sc.grid()
        self.mb_sc.menu = tk.Menu(self.mb_sc)
        self.mb_sc["menu"] = self.mb_sc.menu
        
        #for each spectral clustering algorithm in the list of supported spectral clustering algorithms -> add a button in the menu
        for i, algorithm_clustering in enumerate(self.supported_algorithms_clustering):
            checked_tk = tk.BooleanVar()
            self.mb_sc.menu.add_checkbutton(label = algorithm_clustering, onvalue = True, offvalue = False, variable = checked_tk)
            self.checked.append(checked_tk)            
              
        #number of clusters parameter entry
        ks_insert_label = tk.Label(self.parameters_fr, text="Enter number of clusters for spectral clustering: Enter an int, a list of ints splitted with ',': ")
        ks_insert_label.pack()
        self.ks_tk.pack()
        self.parameters_fr.pack()
        
        #run, back and reset buttons
        run_button = tk.Button(self.sc_fr, text="Run", bg = self.button_bg, fg = self.button_fg, command=self.spectralClusteringRun)
        self.back_button = tk.Button(self.sc_fr, text='Back', bg = self.button_bg, fg = self.button_fg, command = self.back)
        self.reset_button = tk.Button(self.sc_fr, text='Reset Analysis', bg = self.button_bg, fg = self.button_fg, command = self.reset) 

        self.sc_fr.pack()
        self.mb_sc.pack()
            
        run_button.pack()
        self.back_button.pack()
        self.reset_button.pack()    
        
    def communityDetection(self):
        """
        Select the community detection algorithms and, in case of "Asyn FluidC" algorithm select the number of communities to use.
        """
        #if the user press back, return to the "userAnalysisChoice" scene 
        self.previous_scene = "userAnalysisChoice"
        if self.isback or self.isreset:
            #if this page is loaded because the user press the back or the reset button we have to re-initialize some attributes
            self.checked = list()
            self.cd_fr = tk.Frame(self.window)
            self.mb_cd = tk.Menubutton(self.cd_fr, text="Click here and Choose at least one community detection algorithm to use: ", bg = self.button_bg, fg = self.button_fg)
            self.parameters_fr = tk.Frame(self.window)
            self. ks_tk = tk.Entry(self.parameters_fr, text="Number of clusters for clustering", width=100)

        self.analysis_fr.pack_forget()
         
        self.mb_cd.grid()
        self.mb_cd.menu = tk.Menu(self.mb_cd)
        self.mb_cd["menu"] = self.mb_cd.menu
        
        #for each community detection algorithm in the list of supported community detection algorithms -> add a button in the menu
        for i, algorithm_community in enumerate(self.supported_algorithms_communities):
            checked_tk = tk.BooleanVar()
            self.mb_cd.menu.add_checkbutton(label = algorithm_community, onvalue = True, offvalue = False, variable = checked_tk)
            self.checked.append(checked_tk)            
        
        #number of communities parameter entry, used only if Async FluidC algorithm is selected        
        ks_insert_label = tk.Label(self.parameters_fr, text="Enter number of communities for Async FluidC algorithm: Enter an int or a list of ints splitted with ',': ")
        ks_insert_label.pack()
        self.ks_tk.pack()
        self.parameters_fr.pack()
        
        #run, back and reset buttons
        run_button = tk.Button(self.cd_fr, text="Run", bg = self.button_bg, fg = self.button_fg, command = self.communityDetectionRun)
        self.back_button = tk.Button(self.cd_fr, text='Back', bg = self.button_bg, fg = self.button_fg, command = self.back)
        self.reset_button = tk.Button(self.cd_fr, text='Reset Analysis', bg = self.button_bg, fg = self.button_fg, command = self.reset)
         
        self.cd_fr.pack()
        self.mb_cd.pack()
         
        run_button.pack()
        self.back_button.pack()
        self.reset_button.pack()   

    def embeddingClustering(self):
        """
        Select the embedding + clustering algorithms and the parameters to use.
        """
        #if the user press back, return to the "userAnalysisChoice" scene 
        self.previous_scene = "userAnalysisChoice"

        if self.isback or self.isreset:
            #if this page is loaded because the user press the back or the reset button we have to re-initialize some attributes
            self.checked = list()
            self.ec_fr = tk.Frame(self.window)
            self.parameters_fr = tk.Frame(self.window)
            self.ks_tk = tk.Entry(self.parameters_fr, text="Number of clusters for clustering", width=100)
            self.d_tk = tk.Entry(self.parameters_fr, text="d-dimensional embedding: ")
            self.beta_tk = tk.Entry(self.parameters_fr, text="decay factor HOPE embedding: ")
            self.num_walks_tk = tk.Entry(self.parameters_fr, text="Random walk lenght: ")
            self.walk_len_tk = tk.Entry(self.parameters_fr, text="Number of walks per node: ")
            self.mb_ec = tk.Menubutton(self.ec_fr, text="Click here and Choose at least one embedding+clustering algorithm to use: ", bg = self.button_bg, fg = self.button_fg)

        self.analysis_fr.pack_forget()
        
        self.mb_ec.grid()
        self.mb_ec.menu = tk.Menu(self.mb_ec)
        self.mb_ec["menu"] = self.mb_ec.menu
        
        #for each embedding + clustering algorithm in the list of supported embedding + clustering algorithms -> add a button in the menu
        for i, algorithm_emb_clustering in enumerate(self.supported_algorithms_embeddings):
            checked_tk = tk.BooleanVar()
            self.mb_ec.menu.add_checkbutton(label = algorithm_emb_clustering, onvalue = True, offvalue = False, variable = checked_tk)
            self.checked.append(checked_tk)  
        
        #number of clusters parameter entry           
        ks_insert_label = tk.Label(self.parameters_fr, text="Enter number of clusters for embedding + clustering algorithm: Enter an int or a list of ints splitted with ',': ")
        ks_insert_label.pack()
        self.ks_tk.pack()
        
        #d parameter entry for d-dimensional embedding
        d_insert_label = tk.Label(self.parameters_fr, text=" Enter d parameter for d-dimensional embedding: ")
        d_insert_label.pack()
        self.d_tk.pack()
        
        #beta parameter entry for HOPE embedding
        beta_insert_label = tk.Label(self.parameters_fr, text="Enter beta parameter for HOPE embedding: ")
        beta_insert_label.pack()
        self.beta_tk.pack()
        
        #number of walks and walk lenght for node2vec embedding
        num_walks_insert_label = tk.Label(self.parameters_fr, text="Enter the lenght of each random walk: ")
        num_walks_insert_label.pack()
        self.num_walks_tk.pack()
         
        walk_len_insert_label = tk.Label(self.parameters_fr, text="Enter the number of walks per node: ")
        walk_len_insert_label.pack()
        self.walk_len_tk.pack()
        
        self.parameters_fr.pack()
         
        #run, back and reset buttons
        run_button = tk.Button(self.ec_fr, text="Run", bg = self.button_bg, fg = self.button_fg, command = self.embeddingClusteringRun)
        self.back_button = tk.Button(self.ec_fr, text='Back', bg = self.button_bg, fg = self.button_fg, command = self.back)
        self.reset_button = tk.Button(self.ec_fr, text='Reset Analysis', bg = self.button_bg, fg = self.button_fg, command = self.reset)
              
        self.ec_fr.pack()
        self.mb_ec.pack()
        run_button.pack()
        self.back_button.pack()
        self.reset_button.pack()   
            
    def centralityMeasureRun(self):
        """
        For each protein, run all the node centrality algorithms selected.
        """
        cm_touse = []
        
        #check if a supported centrality algorithm is selected
        #checked is a list of booleans, True if the algorithm in position [i] is selected by the user.
        for i, check_button in enumerate(self.checked):
            if check_button.get(): #checked
                cm_touse.append(self.supported_centralities_measures[i])
                
        filepaths = []
        self.cm_fr.pack_forget()
        #for each protein
        for protein, adj in self.adj_matrixs.items():
            
            protein_path = self.proteins_path+protein+".pdb"
            G = nx.from_numpy_matrix(adj)  
            residue_names = self.proteins_residue_names[protein]
            residue_names_1 = np.array(residue_names[:, 1], dtype = str)  
            self.run_fr.pack()
            #for each protein and for each node centrality algorithm to use -> compute, save and plot centralities
            for centrality_choice in cm_touse:
                
                method_to_call = getattr(pcn_final, centrality_choice) 

                compute_tk = tk.Label(self.run_fr, text="Compute {} node centrality for protein {}...".format(centrality_choice, protein)) 
                compute_tk.pack()
                self.window.update()
                centrality_measures = method_to_call(G, residue_names_1)#call the supported method from the pcn_final file
                
                pcn_final.save_centralities(self.output_path, centrality_measures, protein, centrality_choice) #save a txt file 

                plot_tk = tk.Label(self.run_fr, text="Plot {} node centrality for protein {}".format(centrality_choice, protein))
                plot_tk.pack()
                self.window.update()
                pcn_pymol_scripts.pymol_plot_centralities(self.output_path, centrality_measures, protein_path, centrality_choice, self.run_fr, self.window) #plot and save centralities with pymol
                filepath = "{}Centralities{}{}{}Sessions{}{}_{}_session.pse".format(self.output_path, self.add_slash_to_path, centrality_choice, self.add_slash_to_path, self.add_slash_to_path, protein, centrality_choice)
                filepaths.append(filepath)
        
        #next scene
        self.showResults(filepaths)
        
    def spectralClusteringRun(self):
        """
        For each protein and for each number of clusters k, run all the spectral clustering algorithms selected.
        """
        sc_touse = []
                      
        #check if a supported spectral clustering algorithm is selected
        #checked is a list of booleans, True if the algorithm in position [i] is selected by the user.
        for i, check_button in enumerate(self.checked):
            if check_button.get(): #checked
                sc_touse.append(self.supported_algorithms_clustering[i])        
        self.sc_fr.pack_forget()        
        #read k parameter from user entry input
        k_choice = str(self.ks_tk.get())
        if k_choice.split(','):    
            self.ks = [int(k) for k in k_choice.replace(" ","").split(",")] #strip space
      
        filepaths = []
        self.parameters_fr.pack_forget()
        self.run_fr.pack()
        #for each protein
        for protein, adj in self.adj_matrixs.items():
            G = nx.from_numpy_matrix(adj)
            protein_path = self.proteins_path+protein+".pdb"
            residue_names = self.proteins_residue_names[protein]
            residue_names_1 = np.array(residue_names[:, 1], dtype = str)
            #for each spectral clustering algorithm selected
            for algorithm_choice in sc_touse:
                #for each protein, for each spectral clustering algorithm selected and for each number of clusters selected -> compute, save and plot clusters
                for k in self.ks:
                    
                    method_to_call = getattr(pcn_final, algorithm_choice)
                    
                    compute_tk = tk.Label(self.run_fr, text="Compute {} spectral clustering with k = {} for protein {}...".format(algorithm_choice, k, protein)) 
                    compute_tk.pack()
                    self.window.update()
                    labels = method_to_call(adj, n_clusters=k)
                    pcn_final.save_labels(self.output_path, labels, residue_names, protein, algorithm_choice, self.d, self.beta, self.walk_len, self.num_walks)
                    
                    plot_tk = tk.Label(self.run_fr, text="Plot {} spectral clustering with k = {} for protein {}...".format(algorithm_choice, k, protein))
                    plot_tk.pack()
                    self.window.update()
                    pcn_pymol_scripts.pymol_plot(protein_path, self.output_path, "Clusters", algorithm_choice, k, self.run_fr, self.window)
                    filepath = "{}{}{}Sessions{}{}_{}_{}_{}{}_session.pse".format(self.output_path, algorithm_choice, self.add_slash_to_path, self.add_slash_to_path, protein, "Clusters", algorithm_choice, "k", k)
                    filepaths.append(filepath)
                    
                    plot_p_tk = tk.Label(self.run_fr, text="Compute and Plot partecipation coefficients with {} spectral clustering and k = {} for protein {}...".format(algorithm_choice, k, protein))
                    plot_p_tk.pack()
                    self.window.update()
                    p = pcn_final.participation_coefs(G, labels, residue_names_1)
                    pcn_final.save_part_coef(self.output_path, p, protein, algorithm_choice, k)
                    output_path_p = "{}{}{}{}".format(self.output_path, self.add_slash_to_path, algorithm_choice, self.add_slash_to_path)
                    pcn_pymol_scripts.pymol_plot_part_coefs(p, protein_path, output_path_p, algorithm_choice, k, self.run_fr, self.window)
                    filepath_p = "{}Part_coefs_Sessions{}{}_part_coefs_{}_k{}_session.pse".format(output_path_p, self.add_slash_to_path, protein, algorithm_choice, k)
                    filepaths.append(filepath_p)
        
        #next scene
        self.showResults(filepaths)
        
    def communityDetectionRun(self):
        """
        For each protein, run all the community detection algorithms selected with the given parameters.
        """
        cd_touse = []
        #check if a supported community detection algorithm is selected
        #checked is a list of booleans, True if the algorithm in position [i] is selected by the user.
        for i, check_button in enumerate(self.checked):
            if check_button.get(): #checked
                cd_touse.append(self.supported_algorithms_communities[i])
        self.cd_fr.pack_forget()
        
        #if the user wants to use "Asyn FluidC" community detection algorithm, the user have to provide the number of communities to extract
        if "asyn_fluidc" in cd_touse:
            #read k parameter from user entry input
            k_choice = str(self.ks_tk.get())
            if k_choice.split(','):    
                self.ks = [int(k) for k in k_choice.replace(" ","").split(",")] #strip space
        
        self.parameters_fr.pack_forget()
        
        filepaths = []
        #for each protein
        for protein, adj in self.adj_matrixs.items():
            protein_path = self.proteins_path+protein+".pdb"
            G = nx.from_numpy_matrix(adj) 
            residue_names = self.proteins_residue_names[protein]
            residue_names_1 = np.array(residue_names[:, 1], dtype = str)
            
            self.run_fr.pack()
            #for each community detection algorithm
            for algorithm_choice in cd_touse:
                method_to_call = getattr(pcn_final, algorithm_choice)
                       
                #if the algorithm is asyn fluidic
                if (algorithm_choice == 'asyn_fluidc'):
                    #for each protein and for each number of communities k to try -> compute, save and plot communities extracted with Asyn FluidC
                    for k in self.ks:
                                
                        compute_tk = tk.Label(self.run_fr, text="Compute Asyn FluidC communities with k = {} for protein {}...".format(k, protein)) 
                        compute_tk.pack()
                        self.window.update()
                        labels = method_to_call(G, k) #call the method
                        pcn_final.save_labels(self.output_path, labels, residue_names, protein,  method=algorithm_choice) #save the communities as txt file
                        
                        plot_tk = tk.Label(self.run_fr, text="Plot Asyn FluidC with k = {} for protein {}...".format(k, protein))
                        plot_tk.pack()
                        self.window.update()
                        pcn_pymol_scripts.pymol_plot(protein_path, self.output_path, "Communities", algorithm_choice, k, self.run_fr, self.window) #plot and save the communities with pymol
                        filepath = "{}{}{}Sessions{}{}_{}_{}_{}{}_session.pse".format(self.output_path, algorithm_choice, self.add_slash_to_path, self.add_slash_to_path, protein, "Communities", algorithm_choice, "ncoms", k)
                        filepaths.append(filepath)
                        plot_p_tk = tk.Label(self.run_fr, text="Compute and Plot partecipation coefficients with Asyn FluidC and k = {} for protein {}...".format(k, protein))
                        plot_p_tk.pack()
                        self.window.update()                   
                        p = pcn_final.participation_coefs(G, labels, residue_names_1)
                        pcn_final.save_part_coef(self.output_path, p, protein, algorithm_choice, k) #save the part coefs as txt file
                        output_path_p = "{}{}{}".format(self.output_path, algorithm_choice, self.add_slash_to_path)
                        pcn_pymol_scripts.pymol_plot_part_coefs(p, protein_path, output_path_p, algorithm_choice, k, self.run_fr, self.window) #plot and save part coefs with pymol
                        filepath_p = "{}Part_coefs_Sessions{}{}_part_coefs_{}_k{}_session.pse".format(output_path_p, self.add_slash_to_path, protein, algorithm_choice, k)
                        filepaths.append(filepath_p)
                
                else:#if the community detection algorithm is not Asyn Fluidc, no need to specify the number of communities
                    
                    compute_tk = tk.Label(self.run_fr, text="Compute {} communities for protein {}...".format(algorithm_choice, protein)) 
                    compute_tk.pack()
                    self.window.update()
                    labels = method_to_call(G) #call the method 
                    n_coms = int( max(labels) + 1)
                    pcn_final.save_labels(self.output_path, labels, residue_names, protein,  method=algorithm_choice) #save communities as txt 
                    
                    plot_tk = tk.Label(self.run_fr, text="Plot {} with ncoms = {} for protein {}...".format(algorithm_choice, n_coms, protein))
                    plot_tk.pack()
                    self.window.update()
                    pcn_pymol_scripts.pymol_plot(protein_path, self.output_path, "Communities", algorithm_choice, n_coms, self.run_fr, self.window) #plot and save communities with pymol
                    filepath = "{}{}{}Sessions{}{}_{}_{}_{}{}_session.pse".format(self.output_path, algorithm_choice, self.add_slash_to_path, self.add_slash_to_path, protein, "Communities", algorithm_choice, "ncoms", n_coms)
                    filepaths.append(filepath)
                    
                    plot_p_tk = tk.Label(self.run_fr, text="Compute and Plot partecipation coefficients with {} and ncoms = {} for protein {}...".format(algorithm_choice, n_coms, protein))
                    plot_p_tk.pack()
                    self.window.update()                   
                    p = pcn_final.participation_coefs(G, labels, residue_names_1)
                    pcn_final.save_part_coef(self.output_path, p, protein, algorithm_choice, n_coms)
                    output_path_p = "{}{}{}".format(self.output_path, algorithm_choice, self.add_slash_to_path)
                    pcn_pymol_scripts.pymol_plot_part_coefs(p, protein_path, output_path_p, algorithm_choice, n_coms, self.run_fr, self.window)
                    filepath_p = "{}Part_coefs_Sessions{}{}_part_coefs_{}_k{}_session.pse".format(output_path_p, self.add_slash_to_path, protein, algorithm_choice, n_coms)
                    filepaths.append(filepath_p)
        
        #next scene
        self.showResults(filepaths)
        
    def embeddingClusteringRun(self):
        """
        For each protein, run all the embeddings + clustering algorithms selected with the given parameters.
        """    
        ec_touse = []
        #check if a supported embedding + clustering algorithm is selected
        #checked is a list of booleans, True if the algorithm in position [i] is selected by the user.
        for i, check_button in enumerate(self.checked):
            if check_button.get(): #checked
                ec_touse.append(self.supported_algorithms_embeddings[i])
        self.ec_fr.pack_forget()
        
        #read k parameter from user entry input
        k_choice = str(self.ks_tk.get())
        if k_choice.split(','):    
            self.ks = [int(k) for k in k_choice.replace(" ","").split(",")] #strip space
        #read d paramter from user entry input
        self.d = int(self.d_tk.get())
                
        filepaths = []
        self.parameters_fr.pack_forget()
        
        #for each protein 
        for protein, adj in self.adj_matrixs.items():
            protein_path = self.proteins_path+protein+".pdb"
            G = nx.from_numpy_matrix(adj) 
            residue_names = self.proteins_residue_names[protein]
            residue_names_1 = np.array(residue_names[:, 1], dtype = str)
            #for each embedding + clustering algorithm selected
            for algorithm_choice in ec_touse:
                
                #if hope embedding
                if "hope" in algorithm_choice:
                    self.beta = float(self.beta_tk.get())
                #if node2vec embedding
                if "node2vec" in algorithm_choice:
                    self.num_walks = int(self.num_walks_tk.get())
                    self.walk_len =  int(self.walk_len_tk.get())
                
                #for each protein, for each embedding + clustering algorithm selected and for each number of clusters k -> compute, save and plot clusters
                for k in self.ks:

                    method_to_call = getattr(pcn_final, algorithm_choice)
                    compute_tk = tk.Label(self.run_fr, text="Compute {} embedding + clustering with k = {}, d = {},  beta = {}, num_walks = {} and walk_len = {} for protein {}...".format(algorithm_choice, k, self.d, self.beta, self.num_walks, self.walk_len, protein)) 
                    compute_tk.pack()
                    self.run_fr.pack()
                    self.window.update()
                    labels = method_to_call(adj, n_clusters=k, d=self.d, beta=self.beta, walk_len=self.walk_len, num_walks=self.num_walks)
                    pcn_final.save_labels(self.output_path, labels, residue_names, protein, algorithm_choice, self.d, self.beta, self.walk_len, self.num_walks)
                    
                    plot_tk = tk.Label(self.run_fr, text="Plot {} embedding + clustering with k = {}, d = {},  beta = {}, num_walks = {} and walk_len = {} for protein {}...".format(algorithm_choice, k, self.d, self.beta, self.num_walks, self.walk_len, protein)) 
                    plot_tk.pack()
                    self.run_fr.pack()
                    self.window.update()
                    pcn_pymol_scripts.pymol_plot_embeddings(protein_path, self.output_path, "ClustersEmbeddings", algorithm_choice, k, self.d, self.beta, self.walk_len, self.num_walks, self.run_fr, self.window)
                                    
                    if "hope" in algorithm_choice:
                        filepath = "{}{}{}Sessions{}{}_{}_d{}_beta{}_k{}_session.pse".format(self.output_path, algorithm_choice, self.add_slash_to_path, self.add_slash_to_path, protein, algorithm_choice, self.d, self.beta, k)
                        
                    elif "node2vec" in algorithm_choice:
                        filepath = "{}{}{}Sessions{}{}_{}_d{}_wl{}_nw{}_k{}_session.pse".format(self.output_path, algorithm_choice, self.add_slash_to_path, self.add_slash_to_path, protein, algorithm_choice, self.d, self.walk_len, self.num_walks, k)

                    else:
                        filepath = "{}{}{}Sessions{}{}_{}_d{}_k{}_session.pse".format(self.output_path, algorithm_choice, self.add_slash_to_path, self.add_slash_to_path, protein, algorithm_choice, self.d, k)
                        
                    filepaths.append(filepath)
                    
                    plot_p_tk = tk.Label(self.run_fr, text="Compute and Plot partecipation coefficients with {}, k = {}, d = {},  beta = {}, num_walks = {} and walk_len = {} for protein {}...".format(algorithm_choice, k, self.d, self.beta, self.num_walks, self.walk_len, protein)) 
                    plot_p_tk.pack()
                    self.run_fr.pack()
                    self.window.update() 
                    p = pcn_final.participation_coefs(G, labels, residue_names_1)
                    pcn_final.save_part_coef(self.output_path, p, protein, algorithm_choice, k)
                    output_path_p = "{}{}{}{}".format(self.output_path, self.add_slash_to_path, algorithm_choice, self.add_slash_to_path)
                    pcn_pymol_scripts.pymol_plot_part_coefs(p, protein_path, output_path_p, algorithm_choice, k, self.run_fr, self.window)
                    filepath_p = "{}Part_coefs_Sessions{}{}_part_coefs_{}_k{}_session.pse".format(output_path_p, self.add_slash_to_path, protein, algorithm_choice, k)
                    filepaths.append(filepath_p)
        
        #next scene
        self.showResults(filepaths)

    def reset(self):
        """
        Re-initialize all the frames, the list of proteins, the list of selected algorithms and the algorithm parameters selected.
        """
        self.isreset = True
        
        self.comp_adj_fr.pack_forget()
        self.analysis_fr.pack_forget()
        self.cm_fr.pack_forget()
        self.cd_fr.pack_forget()
        self.sc_fr.pack_forget()
        self.ec_fr.pack_forget()
        self.parameters_fr.pack_forget()
        self.run_fr.pack_forget()
        self.results_fr.pack_forget()
        
        self.checked = list()
        self.proteins_list = list()
        self.ks = list()
        self.d = None
        self.beta = None
        self.num_walks = None
        self.walk_len = None
        
        self.select_proteins()
        
    def back(self):
        """
        Go back to the previous scene
        """
        
        self.isback = True 
        if self.previous_scene == "select_proteins":
            self.comp_adj_fr.pack_forget()
            self.analysis_fr.pack_forget()
            self.select_proteins()
        if self.previous_scene == "userAnalysisChoice":
            self.cm_fr.pack_forget()
            self.cd_fr.pack_forget()
            self.sc_fr.pack_forget()
            self.ec_fr.pack_forget()
            self.parameters_fr.pack_forget()
            self.userAnalysisChoice()

    def showResults(self, filepaths):
        """
        Show results.
        For each pyMOL session file created, this fuction creates a button. If the user click this button, it will automatically open the pymol session associated.
        """
        self.run_fr.pack_forget()
        
        #for each pyMOL session file -> create a button
        for filepath in filepaths:
           
           filename = os.path.basename(filepath)
           
           label = tk.Label(self.results_fr, text = filename)
           label.pack()
           button = tk.Button(self.results_fr, text = "Open PyMOL session", bg = self.button_bg, fg = self.button_fg, command = lambda:os.startfile(filepath))
           button.pack()
        
        #only reset button, we are at the end 
        self.reset_button = tk.Button(self.results_fr, text = "Start a new analysis", bg = self.button_bg, fg = self.button_fg, command = self.reset)
        self.results_fr.pack()    
        self.window.update()
        
        self.reset_button.pack(side=tk.BOTTOM)
        self.results_fr.pack()    
        self.window.update()
        
        
    def select_proteins(self):
        """
        Scene to select proteins to analyze, and entry min and max distance threshold for PCN computation.
        """
        if self.isback:
            
            self.choice_fr = tk.Frame(self.window)
            self.comp_adj_fr = tk.Frame(self.window)
            self.proteins_tk = tk.Entry(self.choice_fr, text='PDBs codes:', width=100)
            self.min_tk = tk.Entry(self.choice_fr, text='Min threshold:')
            self.max_tk = tk.Entry(self.choice_fr, text='Max threshold:')
            self.proteins_path_tk = tk.Entry(self.choice_fr, text='PDB Directory:', width=80)
            self.adj_filespath_tk = tk.Entry(self.choice_fr, text='Adj Directory:', width=80)
            self.output_path_tk = tk.Entry(self.choice_fr, text="Output Directory:", width=80)
            
        if self.isreset:    
            
            self.choice_fr = tk.Frame(self.window)
            self.analysis_fr = tk.Frame(self.window)
            self.cm_fr = tk.Frame(self.window)
            self.sc_fr = tk.Frame(self.window)
            self.cd_fr = tk.Frame(self.window)
            self.ec_fr = tk.Frame(self.window)
            self.parameters_fr = tk.Frame(self.window)
            self.comp_adj_fr = tk.Frame(self.window)
            self.run_fr = tk.Frame(self.window)
            self.results_fr = tk.Frame(self.window)
        
            self.proteins_tk = tk.Entry(self.choice_fr, text='PDBs codes:', width=100)
            self.min_tk = tk.Entry(self.choice_fr, text='Min threshold:')
            self.max_tk = tk.Entry(self.choice_fr, text='Max threshold:')

        self.initial_fr_button.pack_forget()
        if(os.path.isfile(os.getcwd()+self.add_slash_to_path+"config.ini")):
              
            config_obj = configparser.ConfigParser()
            config_obj.read(os.getcwd()+self.add_slash_to_path+"config.ini")
            paths = config_obj["user_paths"]
            self.output_path = paths["output_path"]
            self.proteins_path = paths["proteins_path"]
            self.adj_filespath = paths["adj_filespath"]
            paths = tk.Label(self.choice_fr, text="Paths in the config file: \n output path = {} \n proteins_path = {} \n adj_filespath = {}".format(self.output_path, self.proteins_path, self.adj_filespath))
            paths.pack(pady=0)
          
        else:
            
            self.config = False
            
            proteins_path_label = tk.Label(self.choice_fr, text = "Insert Proteins PDB directory:")
            proteins_path_label.pack()
            self.proteins_path_tk.pack()
            
            adj_path_label = tk.Label(self.choice_fr, text = "Insert Adjacency matrixs directory:")
            adj_path_label.pack()
            self.adj_filespath_tk.pack()
            
            output_path_label = tk.Label(self.choice_fr, text = "Insert Output directory:")
            output_path_label.pack()
            self.output_path_tk.pack()

        proteins_insert_label = tk.Label(self.choice_fr, text="Please Insert Protein PDB Identifiers, separated by comma, without .pdb, e.g. 7nxc for 7nxc.pdb")
        proteins_insert_label.pack()
        self.proteins_tk.pack()
          
        adj_pdb_choice = tk.Label(self.choice_fr, text="Format File Input: PDB structures or Preprocessed PCN")
        adj_pdb_choice.pack()
         
        choices = ("pdb", "adj")
        self.choiceVar.set(choices[0])
        cb = ttk.Combobox(self.choice_fr, textvariable=self.choiceVar, values=choices)
        cb.pack()
          
        min_insert_label = tk.Label(self.choice_fr, text="Please enter non covalent bonds threshold distance for PCN costruction")
        min_insert_label.pack()  
        self.min_tk.pack()

        max_insert_label = tk.Label(self.choice_fr, text="Please enter only significant bonds threshold distance for PCN costruction")
        max_insert_label.pack()
        self.max_tk.pack()
          
        submit_button = tk.Button(
                                self.choice_fr,
                                text="Start compute PCNs",
                                bg = self.button_bg, 
                                fg = self.button_fg,
                                command = self.computeOrReadPCN
                                 )
        submit_button.pack(pady=12)
        self.choice_fr.pack()