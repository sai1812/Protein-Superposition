# Protein-Superposition
This repository contains a protein superposition program written in python- Mark42.py which takes two proteins' PDB IDs (in UPPER_case) as inputs and outputs RMSD using quaternion rotation. 

Mark85.py does Multiple pairwise superposition - takes 1 pdb id and 1 text file (with other protein PDB IDs line by line) as inputs and outpts RMSD in form of a sorted table

Note: 
1. Both Mark42.py and Mark85.py uses clustalw for sequence alignemnt. So clustalw should be installed in one's computer before execution.
2. PDB files should be present in the PDB folder and the folder should be present in the specified location (given in the code) 

