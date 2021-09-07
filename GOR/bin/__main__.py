#!/usr/bin/python
#!python


# File Handling . . . . . . . . . . . . . . . . . . . . . . . . . .  . . . . . . . . . . . . . . . . . . . . . .

import sys
from __predictor__ import *

prot_file = open (sys.argv[1], "r")
file_rl = prot_file.readlines()
prot_seq = ""
head = ""
for line in file_rl:
	if line[0:1] == ">":
		head += line.strip()
	elif line[0:1] != ">":
		prot_seq += line.strip()

lis_ps = list(prot_seq)

# Run Predictor . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

helix_propensity = Predictor(helix_table, lis_ps)
beta_propensity = Predictor(extended_table, lis_ps)
coil_propensity = Predictor(coil_table, lis_ps)
turns_propensity = Predictor(turns_table, lis_ps)

# Compare propensity . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

str_data = ""
for i in range(0,len(lis_ps)):
	temp_ls = [helix_propensity[i], beta_propensity[i],coil_propensity[i],turns_propensity[i]]
	temp_str = (temp_ls.index(max(temp_ls)))
	if temp_str == 0:
		str_data += "H"
	elif temp_str == 1:
		str_data += "B" 
	elif temp_str == 2:
		str_data += "C"
	elif temp_str == 3:
		str_data += "T"
# Output . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
print("\nGOR (Garnier–Osguthorpe–Robson) Method\n")
print("Reference: Garnier, J., Osguthorpe, D. J., & Robson, B. (1978).Analysis of\n\
the accuracy and implications of simple methods for predicting the secondary\n\
structure of globular proteins. Journal of molecular biology,120(1), 97–120.\n\
https://doi.org/10.1016/0022-2836(78)90297-8\n")
print("\nProtein: "+head[1:]+"\n")
print("Length: "+str(len(prot_seq))+"\n")

for i in range (0, len(prot_seq),50):
	count = i
	print("\nStructure\t"+str(count+1)+"\t"+str_data[i:i+50]+"\t"+str(count+len(str_data[i:i+50])))
	print("\t\t\t"+"|"*len(str_data[i:i+50])+"\t")
	print("Amino Acids\t"+str(count+1)+"\t"+prot_seq[i:i+50]+"\t"+str(count+len(str_data[i:i+50]))+"\n\n")
	count +=1

print("H: Alpha Helix Conformation")
print("B: Extended Beta-Sheet Conformation")
print("C: Coil Conformation")
print("T: Turns Conformation\n")
