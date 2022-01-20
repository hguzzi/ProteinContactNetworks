from __future__ import print_function
import ast 
import numpy as np
from pymol import cmd 
import sys 
import os 
import pymol

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
    
def get_colors(selection='', quiet=1):
    
    pymol_color_list = []
    
    for tuplepair in pymol.querying.get_color_indices(selection):
        pymol_color_list.append(tuplepair[0])
    
    #pymol_color_list.sort()
    if not int(quiet): print(pymol_color_list)
    pymol_color_list.remove('black')
    pymol_color_list.remove('white')
    return pymol_color_list

cmd.extend('get_colors',get_colors)
cmd.auto_arg[0]['get_colors']=[lambda: cmd.Shortcut(['""','all']), 'selection=', ',']
cmd.auto_arg[1]['get_colors']=[lambda: cmd.Shortcut(['0']), 'quiet=', '']

def pymol_plot(protein_path, output_path, algorithm_type, algorithm_name, k):
    
    cmd.do("delete {}".format("all"))
    cmd.do("load {}".format(protein_path))
    
    protein = os.path.basename(protein_path)
    protein_name = os.path.splitext(protein)[0]
    
    if (algorithm_type == "Communities"):
        ncoms_or_k = "ncoms"
    else:
        ncoms_or_k = "k"
    
    filename = output_path+"{}{}{}{}{}_{}_{}_{}{}.txt".format(algorithm_name, add_slash_to_path, algorithm_type, add_slash_to_path, protein_name, algorithm_type, algorithm_name, ncoms_or_k, k)
    f = open(filename, "r")
    data = f.read()
    dict_node_comms = ast.literal_eval(data)
    f.close()

    colors = get_colors()
    
    for (residue, label) in (dict_node_comms.items()):
        residue_n, residue_chain = residue.split()
        residue_name = residue_n[:3]
        residue_num = residue_n[3:]
        for i in range(k): 
            if(label==i):
                print(residue + " " + colors[i])
                line="color "+colors[i]+", (resi "+ residue_num + " and chain "+ residue_chain + ")"
                cmd.do(line)
    
    if (not os.path.exists("{}{}{}Png".format(output_path, algorithm_name, add_slash_to_path))):
        os.makedirs("{}{}{}Png".format(output_path, algorithm_name, add_slash_to_path))
    if (not os.path.exists("{}{}{}Sessions".format(output_path, algorithm_name, add_slash_to_path))):
        os.makedirs("{}{}{}Sessions".format(output_path, algorithm_name, add_slash_to_path))
        
    cmd.do("capture")
    cmd.do("save {}{}{}Png{}{}_{}_{}_{}{}.png".format(output_path, algorithm_name, add_slash_to_path, add_slash_to_path, protein_name, algorithm_type, algorithm_name, ncoms_or_k, k))          
    cmd.do("save {}{}{}Sessions{}{}_{}_{}_{}{}_session.pse".format(output_path, algorithm_name, add_slash_to_path, add_slash_to_path, protein_name, algorithm_type, algorithm_name, ncoms_or_k, k))
    cmd.do("delete {}".format(protein))

def pymol_plot_chain(protein_path, output_path, algorithm_type, algorithm_name, k, chains_to_delete):
    
    cmd.do("delete {}".format("all"))
    cmd.do("load {}".format(protein_path))
    
    protein = os.path.basename(protein_path)
    protein_name = os.path.splitext(protein)[0]
    
    if (algorithm_type == "Communities"):
        ncoms_or_k = "ncoms"
    else:
        ncoms_or_k = "k"
    
    filename = output_path+"{}{}{}{}{}_{}_{}_{}{}.txt".format(algorithm_name, add_slash_to_path, algorithm_type, add_slash_to_path, protein_name, algorithm_type, algorithm_name, ncoms_or_k, k)
    f = open(filename, "r")
    data = f.read()
    dict_node_comms = ast.literal_eval(data)
    f.close()

    colors = get_colors()
    
    for (residue, label) in (dict_node_comms.items()):
        residue_n, residue_chain = residue.split()
        residue_name = residue_n[:3]
        residue_num = residue_n[3:]
        for i in range(k): 
            if(label==i):
                print(residue + " " + colors[i])
                line="color "+colors[i]+", (resi "+ residue_num + " and chain "+ residue_chain + ")"
                cmd.do(line)
    
    if (not os.path.exists("{}{}{}Png".format(output_path, algorithm_name, add_slash_to_path))):
        os.makedirs("{}{}{}Png".format(output_path, algorithm_name, add_slash_to_path))
    if (not os.path.exists("{}{}{}Sessions".format(output_path, algorithm_name, add_slash_to_path))):
        os.makedirs("{}{}{}Sessions".format(output_path, algorithm_name, add_slash_to_path))
    
    cmd.do("remove resn NAG")
    for chain in chains_to_delete:
        cmd.do("remove chain {}".format(chain))
    
    cmd.do("capture")
    cmd.do("save {}{}{}Png{}{}_{}_{}_{}{}_chain.png".format(output_path, algorithm_name, add_slash_to_path, add_slash_to_path, protein_name, algorithm_type, algorithm_name, ncoms_or_k, k))          
    cmd.do("save {}{}{}Sessions{}{}_{}_{}_{}{}_chain_session.pse".format(output_path, algorithm_name, add_slash_to_path, add_slash_to_path, protein_name, algorithm_type, algorithm_name, ncoms_or_k, k))
    cmd.do("delete {}".format(protein))
    
def pymol_plot_embeddings(protein_path, output_path, algorithm_type, algorithm_name, k, d, beta=None):

    cmd.do("delete {}".format("all"))
    cmd.do("load {}".format(protein_path))
    
    protein = os.path.basename(protein_path)
    protein_name = os.path.splitext(protein)[0]
     
    if (beta is not None):
        filename = output_path+"{}{}{}{}{}_{}_{}_d{}_beta{}_k{}.txt".format(algorithm_name, add_slash_to_path, algorithm_type, add_slash_to_path, protein_name, algorithm_type, algorithm_name, d, beta, k)
    else:
        filename = output_path+"{}{}{}{}{}_{}_{}_d{}_k{}.txt".format(algorithm_name, add_slash_to_path, algorithm_type, add_slash_to_path, protein_name, algorithm_type, algorithm_name, d, k)
    
    f = open(filename, "r")
    data = f.read()
    dict_node_comms = ast.literal_eval(data)
    f.close()
    
    colors = get_colors()
    
    for (residue, label) in (dict_node_comms.items()):
        residue_n, residue_chain = residue.split()
        residue_name = residue_n[:3]
        residue_num = residue_n[3:]
        for i in range(k): 
            if(label==i):
                print(residue + " " + colors[i])
                line="color "+colors[i]+", (resi "+ residue_num + " and chain "+ residue_chain + ")"
                cmd.do(line)
         
    cmd.do("capture")
    
    if (not os.path.exists("{}{}{}Png".format(output_path, algorithm_name, add_slash_to_path))):
        os.makedirs("{}{}{}Png".format(output_path, algorithm_name, add_slash_to_path))
    if (not os.path.exists("{}{}{}Sessions".format(output_path, algorithm_name, add_slash_to_path))):
        os.makedirs("{}{}{}Sessions".format(output_path, algorithm_name, add_slash_to_path))
        
    if (beta is not None):
        cmd.do("save {}{}{}Png{}{}_{}_d{}_beta{}_k{}.png".format(output_path, algorithm_name, add_slash_to_path, add_slash_to_path, protein_name, algorithm_name, d, beta, k))    
        cmd.do("save {}{}{}Sessions{}{}_{}_d{}_beta{}_k{}_session.pse".format(output_path, algorithm_name, add_slash_to_path, add_slash_to_path, protein_name, algorithm_name, d, beta, k))
    else:
        cmd.do("save {}{}{}Png{}{}_{}_d{}_k{}.png".format(output_path, algorithm_name, add_slash_to_path, add_slash_to_path, protein_name, algorithm_name, d, k))   
        cmd.do("save {}{}{}Sessions{}{}_{}_d{}_k{}_session.pse".format(output_path, algorithm_name, add_slash_to_path, add_slash_to_path, protein_name, algorithm_name, d, k))
    
    cmd.do("delete {}".format(protein))
    
def pymol_plot_centralities(centralities, protein_path, output_path, algorithm_name):
    
    cmd.do("delete {}".format("all"))
    cmd.do("load {}".format(protein_path))
    cmd.do("set specular, off")
    protein = os.path.basename(protein_path)
    protein_name = os.path.splitext(protein)[0]
    
    colors = get_colors()
    
    cmd.do("remove hetatm")
    
    for (residue, cent) in (centralities.items()):
        residue_n, residue_chain = residue.split()
        residue_name = residue_n[:3]
        residue_num = residue_n[3:]
        line="alter (resi "+ str(residue_num) + " and chain "+ residue_chain + "), b = "+ str(cent)
        cmd.do(line)
        
    cmd.do("spectrum b, rainbow")
    cmd.do("ramp_new colorbar, none, [{}, {}], rainbow".format(min(centralities.values()), max(centralities.values())))
    
    if (not os.path.exists("{}Centralities{}{}{}Png".format(output_path, add_slash_to_path, algorithm_name, add_slash_to_path))):
        os.makedirs("{}Centralities{}{}{}Png".format(output_path, add_slash_to_path, algorithm_name, add_slash_to_path))
    if (not os.path.exists("{}Centralities{}{}{}Sessions".format(output_path, add_slash_to_path, algorithm_name, add_slash_to_path))):
        os.makedirs("{}Centralities{}{}{}Sessions".format(output_path, add_slash_to_path, algorithm_name, add_slash_to_path))
        
    cmd.do("capture")
    cmd.do("save {}Centralities{}{}{}Png{}{}_{}.png".format(output_path, add_slash_to_path, algorithm_name, add_slash_to_path, add_slash_to_path, protein_name, algorithm_name))          
    cmd.do("save {}Centralities{}{}{}Sessions{}{}_{}_session.pse".format(output_path, add_slash_to_path, algorithm_name, add_slash_to_path, add_slash_to_path, protein_name, algorithm_name))


def pymol_plot_part_coefs(part_coefs, protein_path, output_path):
    
    cmd.do("delete {}".format("all"))
    cmd.do("load {}".format(protein_path))
    cmd.do("set specular, off")
    protein = os.path.basename(protein_path)
    protein_name = os.path.splitext(protein)[0]
    
    colors = get_colors()
    
    cmd.do("remove hetatm")
    
    for (residue, p) in (part_coefs.items()):
        residue_n, residue_chain = residue.split()
        residue_name = residue_n[:3]
        residue_num = residue_n[3:]
        line="alter (resi "+ str(residue_num) + " and chain "+ residue_chain + "), b = "+ str(p)
        cmd.do(line)
        
    cmd.do("spectrum b, rainbow")
    cmd.do("ramp_new colorbar, none, [{}, {}], rainbow".format(min(part_coefs.values()), max(part_coefs.values())))
    
    if (not os.path.exists("{}Part_coefs_Sessions".format(output_path))):
        os.makedirs("{}Part_coefs_Sessions".format(output_path))
  
    cmd.do("save {}Part_coefs_Sessions{}{}_part_coefs_session.pse".format(output_path, add_slash_to_path, protein_name))
