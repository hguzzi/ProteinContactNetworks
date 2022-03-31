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
import gui_constants

##

def compute_adjs():
    
    pdb_list = proteins_list.copy()
    for i, protein in enumerate(pdb_list):
        if (not protein.endswith('.pdb')):
            pdb_list[i] = protein+'.pdb'
    pcn_final.checkIfFilesExists(pdb_list, "pdb", proteins_path)     

def read_adjs():
    pass

def select_proteins():
    
    initial_button.destroy()
    if(os.path.isfile(os.getcwd()+add_slash_to_path+"config.ini")):
        tk.Label(window, text="config file found!")
        config_obj = configparser.ConfigParser()
        config_obj.read(os.getcwd()+add_slash_to_path+"config.ini")
        paths = config_obj["user_paths"]
        output_path = paths["output_path"]
        proteins_path = paths["proteins_path"]
        adj_filespath = paths["adj_filespath"]
    
    else:
        
        proteins_path_tk = tk.Entry(window, text='PDB Directory:')
        proteins_path_tk.pack()
        proteins_path = str(proteins_path_tk.get())
        adj_filespath_tk = tk.Entry(window, text='Adj Directory:')
        adj_filespath_tk.pack()
        adj_filespath = str(adj_filespath_tk.get())
        output_path_tk = tk.Entry(window, text='PDB Directory:')
        output_path_tk.pack()
        output_path = str(output_path_tk.get())
    
    proteins_path = proteins_path+add_slash_to_path
    paths = tk.Label(window, text="Paths in the config file: \n output path = {} \n proteins_path = {} \n adj_filespath = {}".format(output_path, proteins_path, adj_filespath))
    paths.pack()
    
    proteins_insert_label = tk.Label(window, text="Please Insert Protein PDB Identifiers, separated by comma, without .pdb, e.g. 7nxc for 7nxc.pdb")
    proteins_insert_label.pack()
    proteins_tk = tk.Entry(window, text='PDBs codes:')
    proteins_tk.pack()
    proteins = str(proteins_tk.get())
    proteins_list = [protein.casefold() for protein in proteins.replace(" ","").split(",")] #strip space
    adj_pdb_choice = tk.Label(window, text="Format File Input: PDB structures or Preprocessed PCN")
    adj_pdb_choice.pack()
    
    choice_fr = tk.Frame(window)
    choiceVar = tk.StringVar()
    choices = ("pdb", "adj")
    choiceVar.set(choices[0])

    cb = ttk.Combobox(choice_fr, textvariable=choiceVar, values=choices)
    choice_fr.pack()
    cb.pack()
    
    min_insert_label = tk.Label(choice_fr, text="Please enter non covalent bonds threshold distance for PCN costruction")
    min_insert_label.pack()
    min_tk = tk.Entry(choice_fr, text='Min threshold:')
    min_tk.pack()
    min_ = str(min_tk.get())
    max_insert_label = tk.Label(choice_fr, text="Please enter only significant bonds threshold distance for PCN costruction")
    max_insert_label.pack()
    max_tk = tk.Entry(choice_fr, text='Max threshold:')
    max_tk.pack()
    max_ = str(max_tk.get())
    
    if (choiceVar.get() == "pdb"):
        submit_button = tk.Button(
                            choice_fr,
                            text="Start compute PCNs",
                            command = compute_adjs
                                  )
    else:
        submit_button = tk.Button(
                            choice_fr,
                            text="Start compute PCNs",
                            command = read_adjs
                            )

    submit_button.pack()
    
#
supported_algorithms_clustering = ["unnorm_ssc", "norm_ssc", "unnorm_hsc", "norm_hsc", "hsc_shimalik", "ssc_shimalik", "skl_spectral_clustering"]
supported_algorithms_embeddings = [
                                   "unnorm_ssc_hope", "norm_ssc_hope", "unnorm_hsc_hope", "norm_hsc_hope", "hsc_shimalik_hope", "ssc_shimalik_hope",
                                   
                                   "unnorm_ssc_laplacianeigenmaps", "norm_ssc_laplacianeigenmaps", "unnorm_hsc_laplacianeigenmaps", "norm_hsc_laplacianeigenmaps",
                                   "hsc_shimalik_laplacianeigenmaps", "ssc_shimalik_laplacianeigenmaps", 
                                   
                                   "unnorm_ssc_node2vec", "norm_ssc_node2vec", "unnorm_hsc_node2vec", "norm_hsc_node2vec", "hsc_shimalik_node2vec", 
                                   "ssc_shimalik_node2vec",
                                   ]
supported_algorithms_communities = ["louvain", "leiden", "walktrap", "asyn_fluidc", "greedy_modularity", "infomap", "spinglass"]

supported_centralities_measures = ["closeness", "eigenvector_c", "betweenness", "degree_c"]

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

#

window = tk.Tk()
window.title("PCN-Miner 0.9.1-alpha")
window.rowconfigure(0, minsize=800, weight=1)
window.columnconfigure(1, minsize=800, weight=1)
window.geometry('600x400+50+50')
#window.iconbitmap('./assets/pythontutorial.ico') #icon

welcome = tk.Label(window, text="Protein Contact Network Miner 0.9.1-alpha \n Software Available under CC-BY Licence \n Free for Academic Usage ",
                   foreground="black")  
welcome.pack()

initial_fr_button = tk.Frame(window)
initial_button = tk.Button(
                            initial_fr_button,
                            text="Press Enter Button to Continue",
                            command = select_proteins
)
initial_fr_button.pack()
initial_button.pack()

window.mainloop()

