#!/usr/bin/python
#!python


# Predictor_Function

from __gor_table__ import *

def Predictor(Table, lis_ps):

	amino_a_propensity = []	

# Start . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
	pos_start = [0,1,2,3,4,5,6,7]
	for index_s in pos_start:
		temp_index_start = pos_start[:index_s]
		temp_index_start.extend(list(range(index_s,index_s+9)))
		temp_Hval = 0
		for ind in temp_index_start:
			temp_Hval += Table.get(lis_ps[ind])[ind]
		amino_a_propensity.append(temp_Hval)

#  Middle . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	pos_end = []
	for i in range(0,len(lis_ps)):
		temp_seql = lis_ps[i:i+17]
		if len(temp_seql) == 17:
			temp_Hval = 0
			for j in range(0,17):
				temp_Hval += Table.get(temp_seql[j])[j]
			amino_a_propensity.append(temp_Hval)
		elif len(temp_seql) < 17:
			pos_end.append(i)
	
# End . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	for index_e in pos_end:
		temp_index_start = index_e - 8
		if index_e+8 <= pos_end[-1]:
			pass

		else:
			tls = lis_ps[temp_index_start:]
			count = 0
			temp_Hval = 0		
			for ind in tls:
				temp_Hval += Table.get(ind)[count]
				count+=1
			amino_a_propensity.append(temp_Hval)

	return(amino_a_propensity)



