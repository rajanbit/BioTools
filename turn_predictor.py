"""\
Usage:	python turn_predictor.py <fasta_file>

Options:

	<fasta_file>	protein FASTA file
	--help		print help message
"""


import sys

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
	print("Job: Predicting β-Turn Using Chou-Fasman Algorithm\n")

####################################__CHOU- FASMAN ALGORITHM__###################################

	print("|-------------------{ RUNNING CHOU-FASMAN ALGORITHM }-------------------|\n")
	print("Predicting β-Turn Segments...")
	print("Conditions: Pt > 0.000075 and pi > 1.00 and P_Helix < Pi > P_Beta-Sheet")
	print("Where: Pt = f(i)*f(i+1)*f(i+2)*f(i+3) & Pi = the average value of tetrapeptide for P(Turn)")
	new_seq = seq
	for i in range(len(seq)):
		if (i+4) <= len(seq):
			pi = ((turn.get(seq[i])+turn.get(seq[i+1])+turn.get(seq[i+2])+turn.get(seq[i+3]))/4)
			pt = (fi.get(seq[i])*fi1.get(seq[i+1])*fi2.get(seq[i+2])*fi3.get(seq[i+3]))
			p_alpha = ((helix.get(seq[i])+helix.get(seq[i+1])+helix.get(seq[i+2])+helix.get(seq[i+3]))/4)
			p_beta = ((beta.get(seq[i])+beta.get(seq[i+1])+beta.get(seq[i+2])+beta.get(seq[i+3]))/4)
			if pt > 0.000075 and pi > 1.00 and p_alpha < pi > p_beta:
				new_seq = new_seq.replace(seq[i:i+4], "~~~~")
		
	print("Prediction Done/-\n")
	print("|-------------------{ EXITING CHOU-FASMAN ALGORITHM }-------------------|\n")

###################################################################################################
	print("Generating Output in FASTA Format...\n")
	fasta_predicted_helical_seq = ""
	for i in range(0, len(new_seq), 70):
		fasta_predicted_helical_seq += new_seq[i:i+70]+"\n"
	print("Job Done/-\n\n")
	print("	Output: Protein Sequence With Predicted β-Turn Segments")
	print("	--------------------------------------------------------\n\n")
	print(header+fasta_predicted_helical_seq)

...

# python turn_predictor.py <fasta_file>

### Reference:
### Prevelige, P. Jr. and Fasman, G.D., "Chou-Fasman Prediction of the
### Secondary Structure of Proteins," in Prediction of Protein Structure
### and The Priniciples of Protein Conformation (Fasman, G.D., ed.)
### Plenum Press, New York, pp. 391-416 (1989).
