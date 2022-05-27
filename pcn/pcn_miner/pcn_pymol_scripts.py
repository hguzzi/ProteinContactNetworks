from __future__ import print_function
from ast import literal_eval 
import numpy as np
from pymol import cmd
from pymol.querying import get_color_indices 
import os 
from sys import platform

#tk GUI progress bar
import tkinter as tk
from tkinter import ttk

if platform == "linux" or platform == "linux2":
    # linux
    add_slash_to_path = '/'
elif platform == "darwin":
    # OS X
    add_slash_to_path = '/'
elif platform == "win32":
    # Windows...
    add_slash_to_path = '\\' 
    
def get_colors(selection='', quiet=1):
    
    pymol_color_list = []
    
    for tuplepair in get_color_indices(selection):
        pymol_color_list.append(tuplepair[0])
    
    #pymol_color_list.sort()
    if not int(quiet): print(pymol_color_list)
    pymol_color_list.remove('black')
    pymol_color_list.remove('white')
    pymol_color_list.remove('dash')
    return pymol_color_list

cmd.extend('get_colors',get_colors)
cmd.auto_arg[0]['get_colors']=[lambda: cmd.Shortcut(['""','all']), 'selection=', ',']
cmd.auto_arg[1]['get_colors']=[lambda: cmd.Shortcut(['0']), 'quiet=', '']

def pymol_plot(protein_path, output_path, algorithm_type, algorithm_name, k, results_fr = None, window = None):
    
    cmd.do("delete {}".format("all"))
    cmd.do("load {}".format(protein_path))
    
    protein = os.path.basename(protein_path)
    protein_name = os.path.splitext(protein)[0]
    
    if (algorithm_type == "Communities"):
        ncoms_or_k = "ncoms"
        sele_name = "Community"
    else:
        ncoms_or_k = "k"
        sele_name = "Cluster" 
    
    filepath = output_path+"{}{}{}{}{}_{}_{}_{}{}.txt".format(algorithm_name, add_slash_to_path, algorithm_type, add_slash_to_path, protein_name, algorithm_type, algorithm_name, ncoms_or_k, k)
    f = open(filepath, "r")
    data = f.read()
    dict_node_comms = literal_eval(data)
    f.close()

    colors = get_colors()
    cmd.do("remove hetatm")
    
    if results_fr is not None:
        pb = ttk.Progressbar(results_fr, orient="horizontal", mode = "determinate", length = 100)
        pb.pack()
        pb["value"] = 0
        label_tk = tk.Label(results_fr, text = "Current progress {}%".format(pb["value"]))
        label_tk.pack()
        window.update()
    
    n = len(list(dict_node_comms.keys()))
    
    for count, (residue, label) in enumerate((dict_node_comms.items())):
        residue_n, residue_chain = residue.split()
        residue_name = residue_n[:3]
        residue_num = residue_n[3:]
        for i in range(k): 
            if(label==i):
                print(residue + " " + colors[i])
                line="color "+colors[i]+", (resi "+ residue_num + " and chain "+ residue_chain + ")"
                cmd.do(line)
                cmd.do("sele {}, resi {} and chain {}, 1, 0, 1".format("{}{}_{}".format(sele_name, label, colors[i]), residue_num, residue_chain))
        
        if results_fr is not None:
            pb["value"]= round((count/n)*100, 2)
            label_tk['text'] = "Current progress {}%".format(pb["value"])  
            pb.pack()
            label_tk.pack()
            window.update()
    
    if results_fr is not None:
        pb["value"]=round((n/n)*100, 2)
        label_tk['text'] = "Current progress {}%".format(pb["value"])  
        pb.pack()
        label_tk.pack()
        window.update()
    
    if (not os.path.exists("{}{}{}Sessions".format(output_path, algorithm_name, add_slash_to_path))):
        os.makedirs("{}{}{}Sessions".format(output_path, algorithm_name, add_slash_to_path))
        
    cmd.do("save {}{}{}Sessions{}{}_{}_{}_{}{}_session.pse".format(output_path, algorithm_name, add_slash_to_path, add_slash_to_path, protein_name, algorithm_type, algorithm_name, ncoms_or_k, k))
    cmd.do("delete {}".format(protein))
    
def pymol_plot_embeddings(protein_path, output_path, algorithm_type, algorithm_name, k, d, beta=None, walk_len=None, num_walks=None, results_fr = None, window = None):

    cmd.do("delete {}".format("all"))
    cmd.do("load {}".format(protein_path))
    
    protein = os.path.basename(protein_path)
    protein_name = os.path.splitext(protein)[0]
     
    if (beta is not None):
        filepath = output_path+"{}{}{}{}{}_{}_{}_d{}_beta{}_k{}.txt".format(algorithm_name, add_slash_to_path, algorithm_type, add_slash_to_path, protein_name, algorithm_type, algorithm_name, d, beta, k)
    
    elif ((num_walks is not None) and (walk_len is not None)):
        filepath = output_path+"{}{}{}{}{}_{}_{}_d{}_wl{}_nw{}_k{}.txt".format(algorithm_name, add_slash_to_path, algorithm_type, add_slash_to_path, protein_name, algorithm_type, algorithm_name, d, walk_len, num_walks, k)
            
    else:
        filepath = output_path+"{}{}{}{}{}_{}_{}_d{}_k{}.txt".format(algorithm_name, add_slash_to_path, algorithm_type, add_slash_to_path, protein_name, algorithm_type, algorithm_name, d, k)
    
    f = open(filepath, "r")
    data = f.read()
    dict_node_comms = literal_eval(data)
    f.close()
    
    colors = get_colors()
    cmd.do("remove hetatm")
    
    if results_fr is not None:
        pb = ttk.Progressbar(results_fr, orient="horizontal", mode = "determinate", length = 100)
        pb.pack()
        pb["value"] = 0
        label_tk = tk.Label(results_fr, text = "Current progress {}%".format(pb["value"]))
        label_tk.pack()
        window.update()
    
    n = len(list(dict_node_comms.keys()))
    for count, (residue, label) in enumerate(dict_node_comms.items()):
        residue_n, residue_chain = residue.split()
        residue_name = residue_n[:3]
        residue_num = residue_n[3:]
        for i in range(k): 
            if(label==i):
                print(residue + " " + colors[i])
                line="color "+colors[i]+", (resi "+ residue_num + " and chain "+ residue_chain + ")"
                cmd.do(line)
                cmd.do("sele {}, resi {} and chain {}, 1, 0, 1".format("Cluster{}_{}".format(label, colors[i]), residue_num, residue_chain))
        
        if results_fr is not None:
            pb["value"]= round((count/n)*100, 2)
            label_tk['text'] = "Current progress {}%".format(pb["value"])  
            pb.pack()
            label_tk.pack()
            window.update()
    
    if results_fr is not None:
        pb["value"]=round((n/n)*100, 2)
        label_tk['text'] = "Current progress {}%".format(pb["value"])  
        pb.pack()
        label_tk.pack()
        window.update()
    
    if (not os.path.exists("{}{}{}Sessions".format(output_path, algorithm_name, add_slash_to_path))):
        os.makedirs("{}{}{}Sessions".format(output_path, algorithm_name, add_slash_to_path))
        
    if (beta is not None):
        cmd.do("save {}{}{}Sessions{}{}_{}_d{}_beta{}_k{}_session.pse".format(output_path, algorithm_name, add_slash_to_path, add_slash_to_path, protein_name, algorithm_name, d, beta, k))
    elif ((num_walks is not None) and (walk_len is not None)):
        cmd.do("save {}{}{}Sessions{}{}_{}_d{}_wl{}_nw{}_k{}_session.pse".format(output_path, algorithm_name, add_slash_to_path, add_slash_to_path, protein_name, algorithm_name, d, walk_len, num_walks, k))
    else:
        cmd.do("save {}{}{}Sessions{}{}_{}_d{}_k{}_session.pse".format(output_path, algorithm_name, add_slash_to_path, add_slash_to_path, protein_name, algorithm_name, d, k))
    
    cmd.do("delete {}".format(protein))
    
def pymol_plot_centralities(output_path, centralities, protein_path, algorithm_name, results_fr = None, window = None):
    
    cmd.do("delete {}".format("all"))
    cmd.do("load {}".format(protein_path))
    cmd.do("set specular, off")
    protein = os.path.basename(protein_path)
    protein_name = os.path.splitext(protein)[0]
    
    colors = get_colors()
    
    cmd.do("remove hetatm")
    
    if results_fr is not None:
        pb = ttk.Progressbar(results_fr, orient="horizontal", mode = "determinate", length = 100)
        pb.pack()
        pb["value"] = 0
        label = tk.Label(results_fr, text = "Current progress {}%".format(pb["value"]))
        label.pack()
        window.update()
        
    n = len(list(centralities.keys()))
    for count, (residue, cent) in enumerate(centralities.items()):
        residue_n, residue_chain = residue.split(" ")
        residue_name = residue_n[:3]
        residue_num = residue_n[3:]
        line="alter (resi "+ str(residue_num) + " and chain "+ residue_chain + "), b = "+ str(cent)
        cmd.do(line)
    
        if results_fr is not None:
            pb["value"]= round((count/n)*100, 2)
            label['text'] = "Current progress {}%".format(pb["value"])  
            pb.pack()
            label.pack()
            window.update()
    
    if results_fr is not None: 
        pb["value"]=round((n/n)*100, 2)
        label['text'] = "Current progress {}%".format(pb["value"])  
        pb.pack()
        label.pack()
        window.update()
            
    cmd.do("spectrum b, rainbow")
    cmd.do("ramp_new colorbar, none, [{}, {}], rainbow".format(min(centralities.values()), max(centralities.values())))
    
    if (not os.path.exists("{}Centralities{}{}{}Sessions".format(output_path, add_slash_to_path, algorithm_name, add_slash_to_path))):
        os.makedirs("{}Centralities{}{}{}Sessions".format(output_path, add_slash_to_path, algorithm_name, add_slash_to_path))
        
    cmd.do("save {}Centralities{}{}{}Sessions{}{}_{}_session.pse".format(output_path, add_slash_to_path, algorithm_name, add_slash_to_path, add_slash_to_path, protein_name, algorithm_name))


def pymol_plot_part_coefs(part_coefs, protein_path, output_path, algorithm_name, k, results_fr = None, window = None):
    
    cmd.do("delete {}".format("all"))
    cmd.do("load {}".format(protein_path))
    cmd.do("set specular, off")
    protein = os.path.basename(protein_path)
    protein_name = os.path.splitext(protein)[0]
    
    colors = get_colors()
    
    cmd.do("remove hetatm")
    
    if results_fr is not None:
        pb = ttk.Progressbar(results_fr, orient="horizontal", mode = "determinate", length = 100)
        pb.pack()
        pb["value"] = 0
        label = tk.Label(results_fr, text = "Current progress {}%".format(pb["value"]))
        label.pack()
        window.update()
    
    n = len(list(part_coefs.keys()))
    for count, (residue, p) in enumerate(part_coefs.items()):
        
        residue_n, residue_chain = residue.split()
        residue_name = residue_n[:3]
        residue_num = residue_n[3:]
        line="alter (resi "+ str(residue_num) + " and chain "+ residue_chain + "), b = "+ str(p)
        cmd.do(line)
            
        if results_fr is not None:
            pb["value"]= round((count/n)*100, 2)
            label['text'] = "Current progress {}%".format(pb["value"])  
            pb.pack()
            label.pack()
            window.update()
    
    if results_fr is not None:
        pb["value"]=round((n/n)*100, 2)
        label['text'] = "Current progress {}%".format(pb["value"])  
        pb.pack()
        label.pack()
        window.update()
    
    cmd.do("spectrum b, rainbow")
    cmd.do("ramp_new colorbar, none, [{}, {}], rainbow".format(min(part_coefs.values()), max(part_coefs.values())))
    
    if (not os.path.exists("{}Part_coefs_Sessions".format(output_path))):
        os.makedirs("{}Part_coefs_Sessions".format(output_path))
  
    cmd.do("save {}Part_coefs_Sessions{}{}_part_coefs_{}_k{}_session.pse".format(output_path, add_slash_to_path, protein_name, algorithm_name, k))
