"""\
Usage:	python ChouFas_predictor.py <fasta_file>

Options:

	<fasta_file>	protein FASTA file
	--help		print help message
"""


import sys
from time import sleep

...

# Chou-Fasman Amino Acids Propensities Value For Helix
helix = {"A":1.45, "C":0.77,"D":0.98,"E":1.53,"F":1.12,"G":0.53,"H":1.24,"I":1.00,"K":1.07,
"L":1.34,"M":1.20,"N":0.73,"P":0.59,"Q":1.17,"R":0.79,"S":0.79,"T":0.82,"V":1.14,"W":1.14,"Y":0.61}

# Chou-Fasman Amino Acids Propensities Value For B-Sheet
beta = {"A":0.97, "C":1.30,"D":0.80,"E":0.26,"F":1.28,"G":0.81,"H":0.71,"I":1.60,"K":0.74,
"L":1.22,"M":1.67,"N":0.65,"P":0.62,"Q":1.23,"R":0.90,"S":0.72,"T":1.20,"V":1.65,"W":1.19,"Y":1.29}

# Chou-Fasman Amino Acids Propensities Value For B-Turn
turn = {"A":0.66, "C":1.19,"D":1.46,"E":0.74,"F":0.60,"G":1.56,"H":0.95,"I":0.47,"K":1.01,
"L":0.59,"M":0.60,"N":1.56,"P":1.52,"Q":0.98,"R":0.95,"S":1.43,"T":0.96,"V":0.50,"W":0.96,"Y":1.14}

fi = {"A":0.060, "C":0.149,"D":0.147,"E":0.056,"F":0.059,"G":0.102,"H":0.140,"I":0.043,"K":0.055,
"L":0.061,"M":0.068,"N":0.161,"P":0.102,"Q":0.074,"R":0.070,"S":0.120,"T":0.086,"V":0.062,"W":0.077,"Y":0.082}

fi1 = {"A":0.076, "C":0.053,"D":0.110,"E":0.060,"F":0.041,"G":0.085,"H":0.047,"I":0.034,"K":0.115,
"L":0.025,"M":0.082,"N":0.083,"P":0.301,"Q":0.098,"R":0.106,"S":0.139,"T":0.108,"V":0.048,"W":0.013,"Y":0.065}

fi2 = {"A":0.035, "C":0.117,"D":0.179,"E":0.077,"F":0.065,"G":0.190,"H":0.093,"I":0.013,"K":0.072,
"L":0.036,"M":0.014,"N":0.191,"P":0.034,"Q":0.037,"R":0.099,"S":0.125,"T":0.065,"V":0.028,"W":0.064,"Y":0.114}

fi3 = {"A":0.058, "C":0.128,"D":0.081,"E":0.064,"F":0.065,"G":0.152,"H":0.054,"I":0.056,"K":0.095,
"L":0.070,"M":0.055,"N":0.091,"P":0.068,"Q":0.098,"R":0.085,"S":0.106,"T":0.079,"V":0.053,"W":0.167,"Y":0.125}

# Function Progress Bar
def progress(message):
	for i in range(51):
		sys.stdout.write('\r')
		sys.stdout.write(message + "[%-50s] %d%%" % ('#'*i, 2*i))
		sys.stdout.flush()
		sleep(0.03)

# Function for Beta-Turn Prediction
def turn_predictor(seq):
	new_seq = seq
	for i in range(len(seq)):
		if (i+4) <= len(seq):
			pi = ((turn.get(seq[i])+turn.get(seq[i+1])+turn.get(seq[i+2])+turn.get(seq[i+3]))/4)
			pt = (fi.get(seq[i])*fi1.get(seq[i+1])*fi2.get(seq[i+2])*fi3.get(seq[i+3]))
			p_alpha = ((helix.get(seq[i])+helix.get(seq[i+1])+helix.get(seq[i+2])+helix.get(seq[i+3]))/4)
			p_beta = ((beta.get(seq[i])+beta.get(seq[i+1])+beta.get(seq[i+2])+beta.get(seq[i+3]))/4)
			if pt > 0.000075 and pi > 1.00 and p_alpha < pi > p_beta:
				new_seq = new_seq.replace(seq[i:i+4], "~~~~")
	for x in new_seq:
		if x != "~":
			new_seq = new_seq.replace(x, "-")
	new_seq = new_seq.replace("~", "T")
	return(new_seq)

# Function for Helix Prediction
def helix_predictor(seq):
	probable_helix_index=[]
	temp_seq1 = ""
	prot_helix_prob = []
	prot_beta_prob = []
	for i in range(0,len(seq)):
		prot_helix_prob.append(helix.get(seq[i]))
		prot_beta_prob.append(beta.get(seq[i]))
		if helix.get(seq[i]) > 1:
			temp_seq1 += "H"
		elif helix.get(seq[i]) < 1:
			temp_seq1 += "-"
	nucl_index = [i for i in range(len(temp_seq1)) if temp_seq1.startswith("HHHH", i)]
	break_index = [i for i in range(len(temp_seq1)) if temp_seq1.startswith("----", i)]
	for nucl in nucl_index:
		n1 = nucl
		n6 = nucl+6
		for brk in break_index:
			if brk > nucl and brk > n6:
				sum_in_index = brk - n6
				if n1-sum_in_index < 0:
					probable_helix_index.append([0,n6+sum_in_index])
				else:
					probable_helix_index.append([n1-sum_in_index,n6+sum_in_index])
			elif brk < nucl:
				sum_in_index = n1-brk-1
				if n6+sum_in_index > len(seq):
					probable_helix_index.append([n1-sum_in_index,len(seq)])
				else:
					probable_helix_index.append([n1-sum_in_index,n6+sum_in_index])
	if break_index == []:
		probable_helix_index.append([0,len(seq)])
	probable_helix_index_drm = [] 
	for x in probable_helix_index:
		if x not in probable_helix_index_drm:
			probable_helix_index_drm.append(x)
	probable_helix_index_f1 = []
	for data in probable_helix_index_drm:
		if "----" not in temp_seq1[data[0]:data[1]]:
			probable_helix_index_f1.append(data)
	for data in probable_helix_index_f1:
		p_alpha = sum(prot_helix_prob[data[0]:data[1]]) / len(prot_helix_prob[data[0]:data[1]])
		p_beta = sum(prot_beta_prob[data[0]:data[1]]) / len(prot_beta_prob[data[0]:data[1]])
		if p_alpha > 1.03 and p_alpha > p_beta:
			seq = seq.replace(seq[data[0]:data[1]],len(seq[data[0]:data[1]])*"~")
	for x in seq:
		if x != "~":
			seq = seq.replace(x, "-")
	seq = seq.replace("~", "H")
	return(seq)

# Function for B-sheet Prediction
def sheet_predictor(seq):
	probable_beta_index=[]
	temp_seq1 = ""
	prot_helix_prob = []
	prot_beta_prob = []
	for i in range(0,len(seq)):
		prot_helix_prob.append(helix.get(seq[i]))
		prot_beta_prob.append(beta.get(seq[i]))
		if beta.get(seq[i]) > 1:
			temp_seq1 += "B"
		elif beta.get(seq[i]) < 1:
			temp_seq1 += "-"
	nucl_index = [i for i in range(len(temp_seq1)) if temp_seq1.startswith("BBB", i)]
	break_index = [i for i in range(len(temp_seq1)) if temp_seq1.startswith("----", i)]
	for nucl in nucl_index:
		n1 = nucl
		n3 = nucl+3
		n5 = nucl+5
		for brk in break_index:
			if brk > nucl and brk > n3:
				sum_in_index = brk - n3
				if n1-sum_in_index < 0:
					probable_beta_index.append([0,n3+sum_in_index])
				else:
					probable_beta_index.append([n1-sum_in_index,n3+sum_in_index])
			elif brk < nucl:
				sum_in_index = n1-brk-1
				if n3+sum_in_index > len(seq):
					probable_beta_index.append([n1-sum_in_index,len(seq)])
				else:
					probable_beta_index.append([n1-sum_in_index,n3+sum_in_index])
			elif brk > nucl and brk > n5:
				sum_in_index = brk - n5
				if n1-sum_in_index < 0:
					probable_beta_index.append([0,n5+sum_in_index])
				else:
					probable_beta_index.append([n1-sum_in_index,n5+sum_in_index])
			elif brk < nucl:
				sum_in_index = n1-brk-1
				if n5+sum_in_index > len(seq):
					probable_beta_index.append([n1-sum_in_index,len(seq)])
				else:
					probable_beta_index.append([n1-sum_in_index,n5+sum_in_index])
	if break_index == []:
		probable_beta_index.append([0,len(seq)])
	probable_beta_index_drm = [] 
	for x in probable_beta_index:
		if x not in probable_beta_index_drm:
			probable_beta_index_drm.append(x)
	probable_beta_index_f1 = []
	for data in probable_beta_index_drm:
		if "----" not in temp_seq1[data[0]:data[1]]:
			probable_beta_index_f1.append(data)
	for data in probable_beta_index_f1:
		p_alpha = sum(prot_helix_prob[data[0]:data[1]]) / len(prot_helix_prob[data[0]:data[1]])
		p_beta = sum(prot_beta_prob[data[0]:data[1]]) / len(prot_beta_prob[data[0]:data[1]])
		if p_beta > 1.05 and p_beta > p_alpha:
			seq = seq.replace(seq[data[0]:data[1]],len(seq[data[0]:data[1]])*"~")
	for x in seq:
		if x != "~":
			seq = seq.replace(x, "-")
	seq = seq.replace("~", "B")
	return(seq)

# Function to remove Overlaps
def grouper(olist):
    prev = None
    group = []
    for item in olist:
        if not prev or item - prev <= 15:
            group.append(item)
        else:
            yield group
            group = [item]
        prev = item
    if group:
        yield group

def overlap(seq1, seq2, seq3):
	seq1 = list(seq1)
	seq2 = list(seq2)
	seq3 = list(seq3)
	nseq = ""
	overlist = []
	ilist = []
	for x in range(0,len(seq1)):
		if seq1[x] == "H" and seq2[x] == "B":
			overlist.append(x)
			nseq += "-"
		elif seq1[x] == "H" and seq2[x] =="-":
			nseq += "H"
		elif seq1[x] == "-" and seq2[x] =="B":
			nseq += "B"
		elif seq1[x] == "-" and  seq2[x] =="-":
			nseq +="-"
	nseq = list(nseq)
	nseq1 = []
	for i in range(0,len(seq3)):
		if seq3[i] == "T":
			nseq1.append("T")
		else:
			nseq1.append(nseq[i])
	nseq2 = "".join(nseq1)
	a = (list(grouper(overlist)))
	overlap_index = []
	for data in a:
		if len(data) > 1:
			overlap_index.append([data[0],data[-1]])
		else:
			overlap_index.append([data[0]])
	prot_helix_prob = []
	prot_beta_prob = []
	for i in range(0,len(seq)):
		prot_helix_prob.append(helix.get(seq[i]))
		prot_beta_prob.append(beta.get(seq[i]))
	for data in overlap_index:
		if len(data) > 1:
			p_alpha = sum(prot_helix_prob[data[0]:data[1]]) / len(prot_helix_prob[data[0]:data[1]])
			p_beta = sum(prot_beta_prob[data[0]:data[1]]) / len(prot_beta_prob[data[0]:data[1]])
			if p_alpha > p_beta:
				nseq2 = nseq2.replace(nseq2[data[0]:data[1]],len(nseq2[data[0]:data[1]])*"H")
			elif p_beta > p_alpha:
				nseq2 = nseq2.replace(nseq2[data[0]:data[1]],len(nseq2[data[0]:data[1]])*"B")
		else:
			pass
	nseq3 = []
	for i in range(0,len(seq3)):
		if seq3[i] == "T":
			nseq3.append("T")
		else:
			nseq3.append(nseq2[i])
	nseq4 = "".join(nseq3)
	return(nseq4)

# argv[1] is not defined
if len(sys.argv) == 1:
	print (__doc__)

# argv[1] is --help
elif sys.argv[1] == "--help":
	print (__doc__)

# argv[1] is a fasta file
else:
# reading the fasta file
	print("\nParsing Protein Sequence...\n")
	file_f = sys.argv[1]
	fasta_file = open(file_f)
	fasta_rec = fasta_file.readlines()

# separating header from sequence
	header =""
	seq =""
	for line in fasta_rec:
		if line[0:1] == ">":
			header = line
		else:
			seq += line.strip()
	print("Sequence Information")
	print("Sequence Origin: "+header.strip()[1:])
	print("Sequence Length:",len(seq),"\n")
	print("Job: Predicting Protein Secondary Structure Using Chou-Fasman Algorithm\n")

# Running Everything
	print("|-------------------{ RUNNING CHOU-FASMAN ALGORITHM }-------------------|\n")

	print("Predicting Helical Segments...")
	print("Condition: { P_alpha > 1.03 and P_alpha > P_beta}\n")
	progress("Helix Prediction | ")
	h = helix_predictor(seq)
	print("\nHelix Prediction | Done\n")

	print("Predicting B-sheet Segments...")
	print("Condition: { P_beta > 1.05 and P_beta > P_alpha}\n")
	progress("B-sheet Prediction | ")
	s = sheet_predictor(seq)
	print("\nB-sheet Prediction | Done\n")

	print("Predicting β-Turn Segments...")
	print("Conditions: Pt > 0.000075 and pi > 1.00 and P_Helix < Pi > P_Beta-Sheet")
	print("Where: Pt = f(i)*f(i+1)*f(i+2)*f(i+3) & Pi = the average value of tetrapeptide for P(Turn)\n")
	progress("B-turn Prediction | ")
	t = turn_predictor(seq)
	print("\nB-turn Prediction | Done\n")

	print("|-------------------{ EXITING CHOU-FASMAN ALGORITHM }-------------------|\n")

	print("Output: Protein Sequence With Predicted Alpha-Helical Segments, Beta-Sheet Segments and Beta-Turns...\n")

	try:
		fseq = overlap(h,s,t)
		fseq = fseq.replace("H",('\033[31m'+"H"+'\033[0m'))
		fseq = fseq.replace("B",('\033[34m'+"B"+'\033[0m'))
		fseq = fseq.replace("T",('\033[33m'+"T"+'\033[0m'))
		fseq = fseq.replace("-",('\033[32m'+"-"+'\033[0m'))
		print(fseq+"\n")
	except:
		fseq = overlap(h,s,t)
		print(fseq+"\n")
...

# python ChouFas_predictor.py <fasta_file>

### Reference:
### Prevelige, P. Jr. and Fasman, G.D., "Chou-Fasman Prediction of the
### Secondary Structure of Proteins," in Prediction of Protein Structure
### and The Priniciples of Protein Conformation (Fasman, G.D., ed.)
### Plenum Press, New York, pp. 391-416 (1989).
