#!/usr/bin/python
#!python

print("""
######################################################################
#                                                                    #
#  Nussinov algorithm developed by Ruth Nussinov is a nucleic acid   #
#  structure prediction algorithm used to predict the secondary      #
#  structure of RNA by base pair maximization.                       #
#                                                                    #
#  The Program construct a Nussinov Matrix from the RNA sequence     #
#  using NumPy library and the folds are predicted by Traceback      #
#  from the upper right corner of the Nussinov Matrix.               #
#                                                                    #
# ------------------------------------------------------------------ #
#                                                                    #
#  REFERENCE                                                         #
#                                                                    #
#  Nussinov, R., & Jacobson, A. B. (1980). Fast algorithm for        #
#  predicting the secondary structure of single-stranded RNA.        #
#  Proceedings of the National Academy of Sciences of the United     #
#  States of America, 77(11), 6309–6313.                             #
#  https://doi.org/10.1073/pnas.77.11.6309                           #
#                                                                    #
######################################################################
""")

# IMPORTING_MODULES

import sys
import numpy as np

# INPUT Sequence

seq = ""
fl = open(sys.argv[1], "r")
data = fl.readlines()
for line in data:
	if line[0:1] == ">":
		head = line
	else:
		seq += line.strip()

# RNA Sequence

rna_seq = seq.replace("T","U")

# Nussinov Matrix Construction

n_matrix = np.full((len(rna_seq), len(rna_seq)), 0)
for i in range(0, len(rna_seq)):
	for j in range(i+1, len(rna_seq)):
		if (rna_seq[i] == "A" and rna_seq[j] =="U") or (rna_seq[i] == "U" and rna_seq[j] =="A"):
			n_matrix[i,j] = 1

		elif (rna_seq[i] == "G" and rna_seq[j] =="C") or (rna_seq[i] == "C" and rna_seq[j] =="G"):
			n_matrix[i,j] = 1

		elif (rna_seq[i] == "G" and rna_seq[j] =="U") or (rna_seq[i] == "U" and rna_seq[j] =="G"):
			n_matrix[i,j] = 1

print(n_matrix)
