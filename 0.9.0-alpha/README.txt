RUN setup.bat: create a clean anaconda enviroment and install all the libraries

FOR WINDOWS USERS: MAKE SURE YOU HAVE INSTALLED git

ISTRUCTION FOR THE RUN OF pcn_main

1.OPEN Anaconda prompt 

2.conda activate PCN 

3.python pcn_main.py

Check if the config.ini file exists. If not exists, create it.

Check if the user want to use the paths saved in the config.ini file. If he wants to change the paths:
	
	Enter the directory where you want to save your Outputs.

	Enter the directory where you have saved the pdb files.
	
	Choose to start with pdb files or with adj matrixs.
	(If you have choose to start with adj matrixs, Enter the directory where are saved the matrixs)

Enter the names of the proteins you want to study, the pdb file of a protein must be saved
in the given Protein Directory (TO ADD: if the pdb of a protein isn't in the Protein Directory, automatically fetch them from the PDB database)

Check if the user want to compute a centrality measure of the PCN. If he wants: compute the centrality measure of the PCN and plot it with pyMOL. 

Choose the type of study you want to do: Spectral Clustering, Community Extraction, Embeddings + Clustering
(if spectral clustering is involved, provide the numbers of clusters ks you want to use)
(if embedding, provide the parameters)

Choose one, a list, or all supported algorithms.

Plot the clusters/communities with pyMOL.