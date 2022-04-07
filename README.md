# Protein Contact Networks Miner: A tool for the Analysis of Protein Contact Networks

Protein Contact Networks Miner is a command line tool designed for annotate allosteric domains of a protein based of his rappresentation trough a unweighted graph of amino acids significant and non covalent interactions, this graph is also called Protein Contact Network. PCN Miner now has a Graphic User Interface. 

Protein Contact Network is an unweighted graph: the nodes of the graphs are the amino acids and exists an edge that connect two nodes (amino acids) i and j only if the euclidean distance between this two amino acids is between 4 Angstrom (only non covalent interactions) and 8 Angstrom (only significant interactions). The distance between two aminoacids i and j is approssimated by the distance between the Alpha Carbon of the amino acids. The user can modify the only covalent (min) and the only significant (max) threshold distance for PCN construction. 

![image](https://user-images.githubusercontent.com/87126937/162151753-43c6157b-028a-45e2-9aeb-dafd912d4162.png)

![image](https://user-images.githubusercontent.com/87126937/162151714-bf5ce554-14ad-4100-b4e9-6d95af19bca0.png)

PCN global (like graph diameter) or local descriptors (like node centrality measures) are useful to model and analyse protein functions. PCN Miner allow to identify modules (also called communities or clusters) in protein molecules using three different approaches: 
  1. spectral clustering: extract clusters from a graph with a clustering approach based on the Laplacian matrix eigenvectors following the guidelines given    in: A tutorial on spectral clustering [1];
  2. embedding+clustering: uses one of the embedding algorithm in the GEM library [2] and then apply spectral clustering;
  3. community detection: uses one of the community detection algorithm in the cdlib library [3].

Supported Algorithms:
  
  1. Spectral Clustering: Both Hard (K-Means) and Soft (Fuzzy C-Means) clustering approach used on the eigenvectors of the Laplacian matrix (both normalized or unnormalized form). Is also supported the Shi Malik spectral clustering approach that resolves the generalized eigenvalues problem;
  2. Embedding + Clustering: Node2vec, HOPE, Laplacianeigenmaps embedding followed by a supported spectral clustering algorithm;
  3. Community Detection:  Louvain, Leiden, Walktrap, Infomap, Asyn FluidC, Greedy Modularity, Spinglass;
  4. Centrality Measures: Closeness Centrality, Betweenness Centrality, Eigenvector Centrality, Degree Centrality.

Outputs (node centrality, clusters or communities) of the supported algorithms are then plotted on the 3D structure of the protein using PyMol scripts.

Third part softwares needed:
  
    -Git: https://git-scm.com/downloads
    -Anaconda3: https://www.anaconda.com/products/individual
    -PyMOL: https://pymol.org/2/#download

Required libraries: numpy, networkx, regex, scipy, fuzzy-c-means, cdlib, GEM, node2vec, pymol, pytz, python-dateutil, pyinstaller.
  
This libraries are automatically installed when the user runs setupWindows.bat or setupLinux-MACOSX.sh

How to install it:

-S.O. Windows:

              git clone https://github.com/hguzzi/ProteinContactNetworks.git
                                      cd 0.9.1-alpha
                                     setupWindows.bat
        
-S.O. Linux-MACOSX:

              git clone https://github.com/hguzzi/ProteinContactNetworks.git
                                      cd 0.9.1-alpha
                                source setupLinux-MacOSX.sh  
    
How to use the command line version:

                                    conda activate PCN
                                      cd 0.9.1-alpha
                                    python pcn_main.py

How to use the GUI version: #TO DO: use pyinstaller to build an .exe file from the pcn_gui.py script
                                    
                                    conda activate PCN
                                      cd 0.9.1-alpha       
                                  python pcn_gui_main.py

Example:
  
Entry PDB code: 6VXX
Description: SARS CoV 2 Spike protein closed form
                                    
Method: Community Detection
Algorithm: Leiden
Number of communities extracted: 20 

![image](https://user-images.githubusercontent.com/87126937/162151095-3ddc1177-3b32-4407-b6d7-06eb4dab9b3e.png)

Method: Centrality Analysis
Algorithm: Eigenvector Centrality

![image](https://user-images.githubusercontent.com/87126937/162151265-a64b2af6-bb15-41eb-883f-a4cc1779439d.png)


References:
  
  [1] von Luxburg, U. A tutorial on spectral clustering. Stat Comput 17, 395–416 (2007). https://doi.org/10.1007/s11222-007-9033-z;
  
  [2] https://github.com/palash1992/GEM;
  
  [3] G. Rossetti, L. Milli, R. Cazabet. CDlib: a Python Library to Extract, Compare and Evaluate Communities from Complex Networks. Applied Network Science Journal. 2019. DOI:10.1007/s41109-019-0165-9 https://github.com/GiulioRossetti/cdlib
