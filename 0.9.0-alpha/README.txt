RUN setup.bat: create a clean anaconda enviroment 
RUN install_libraries.bat: install all the libraries

FOR WINDOWS USERS: MAKE SURE YOU HAVE INSTALLED git

ISTRUCTION FOR THE RUN OF pcn_main

1.OPEN Anaconda prompt 

2.conda activate PCN 

3.python pcn_main.py

Enter the directory where you want to save your Outputs.

Enter the directory where you have saved the pdb files.

Choose to start with pdb files or with adj matrixs.
(If you have choose to start with adj matrixs, Enter the directory where are saved the matrixs)

Enter the names of the proteins you want to study, the pdb file of a protein must be saved
in the given Protein Directory (TO ADD: if the pdb of a protein isn't in the Protein Directory, automatically fetch them from the PDB database)

Choose the type of study you want to do: Spectral Clustering, Community Extraction, Embeddings + Spectral Clustering
(if spectral clustering is involved, provide the numbers of clusters ks you want to use)
(if embedding, provide the parameters)

Choose one of the supported algorithms.