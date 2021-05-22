### importing modules
####################################################################################

import sys
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.text import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties

### Functions used
####################################################################################

def aa_info(aa, x, y, yscale=1, ax=None):
    text = aas[aa]

    t = mpl.transforms.Affine2D().scale(1*globalscale, yscale*globalscale) + \
        mpl.transforms.Affine2D().translate(x,y) + ax.transData
    p = PathPatch(text, lw=0, fc=aa_clr[aa],  transform=t)
    if ax != None:
        ax.add_artist(p)
    return p

def log2(p):
	if p == 0:
		log_aa = 0
	else:
		log_aa = math.log2(p)
	return(log_aa)

file_f = sys.argv[1]
file_r = open(file_f, "r")
fast_f = file_r.readlines()
seq_list = []
for line in fast_f:
	if line[0] != ">":
		seq_list.append(line.strip())
def score_order(score):  
	lst = len(score)  
	for i in range(0, lst):  
		for j in range(0, lst-i-1):  
			if (score[j][1] > score[j + 1][1]):  
				temp = score[j]  
				score[j]= score[j + 1]  
				score[j + 1]= temp  
	return score

### Fasta file handling and bit_score calculation
####################################################################################
rm_ls = ['','\n','\t']
for item in rm_ls:
	while item in seq_list: 
		seq_list.remove(item)
aa_scores = []
for i in range(0,len(seq_list[0])):
	A = C = D = E = F = G = H = I = K = L = M = N = P = Q = R = S = T = W = Y = V = O = U = 0
	for data in seq_list:
		if data[i] == "A":
			A += 1
		elif data[i] == "C":
			C += 1
		elif data[i] == "D":
			D += 1
		elif data[i] == "E":
			E += 1
		elif data[i] == "F":
			F += 1
		elif data[i] == "G":
			G += 1
		elif data[i] == "H":
			H += 1
		elif data[i] == "I":
			I += 1
		elif data[i] == "K":
			K += 1
		elif data[i] == "L":
			L += 1
		elif data[i] == "M":
			M += 1
		elif data[i] == "N":
			N += 1
		elif data[i] == "P":
			P += 1
		elif data[i] == "Q":
			Q += 1
		elif data[i] == "R":
			R += 1
		elif data[i] == "S":
			S += 1
		elif data[i] == "T":
			T += 1
		elif data[i] == "W":
			W += 1
		elif data[i] == "Y":
			Y += 1
		elif data[i] == "V":
			V += 1
		elif data[i] == "O":
			O += 1
		elif data[i] == "U":
			U += 1

	tot = A + C + D + E + F + G + H + I + K + L + M + N + P + Q + R + S + T + W + Y + V + O + U

	pA = A/tot
	pC = C/tot
	pD = D/tot
	pE = E/tot
	pF = F/tot
	pG = G/tot
	pH = H/tot
	pI = I/tot
	pK = K/tot
	pL = L/tot
	pM = M/tot
	pN = N/tot
	pP = P/tot
	pQ = Q/tot
	pR = R/tot
	pS = S/tot
	pT = T/tot
	pW = W/tot
	pY = Y/tot
	pV = V/tot
	pO = O/tot
	pU = U/tot

	i = (math.log2(4.0)) + (pA*log2(pA)) + (pC*log2(pC)) +(pD*log2(pD)) + (pE*log2(pE)) + (pF*log2(pF)) + (pG*log2(pG)) +(pH*log2(pH)) + (pI*log2(pI)) +(pK*log2(pK)) + (pL*log2(pL)) +(pM*log2(pM)) + (pN*log2(pN)) + (pP*log2(pP)) + (pQ*log2(pQ)) +(pR*log2(pR)) + (pS*log2(pS)) + (pT*log2(pT)) + (pW*log2(pW)) +(pY*log2(pY)) + (pV*log2(pV)) +(pO*log2(pO)) + (pU*log2(pU))

	score = [("A", pA*i), ("C", pC*i), ("D", pD*i), ("E", pE*i), ("F", pF*i), ("G", pG*i), ("H", pH*i), ("I", pI*i), ("K", pK*i), ("L", pL*i), ("M", pM*i), ("N", pN*i), ("P", pP*i), ("Q", pQ*i), ("R", pR*i), ("S", pS*i), ("T", pT*i), ("W", pW*i), ("Y", pY*i), ("V", pV*i),("O", pO*i), ("U", pU*i)]

	aa_scores.append(score_order(score))

	A = C = D = E = F = G = H = I = K = L = M = N = P = Q = R = S = T = W = Y = V = O = U = 0

### plotting the bit_score for amino acids
####################################################################################

font_p = FontProperties(weight="bold")

globalscale = 1.35

aas = {"G" : TextPath((-0.3, 0), "G", size=1, prop=font_p),
            "S" : TextPath((-0.3, 0), "S", size=1, prop=font_p),
            "T" : TextPath((-0.3, 0), "T", size=1, prop=font_p),
            "Y" : TextPath((-0.3, 0), "Y", size=1, prop=font_p),
            "C" : TextPath((-0.3, 0), "C", size=1, prop=font_p),
            "Q" : TextPath((-0.3, 0), "Q", size=1, prop=font_p),
            "N" : TextPath((-0.3, 0), "N", size=1, prop=font_p),
            "K" : TextPath((-0.3, 0), "K", size=1, prop=font_p),
            "R" : TextPath((-0.3, 0), "R", size=1, prop=font_p),
            "H" : TextPath((-0.3, 0), "H", size=1, prop=font_p),
            "D" : TextPath((-0.3, 0), "D", size=1, prop=font_p),
            "E" : TextPath((-0.3, 0), "E", size=1, prop=font_p),
            "A" : TextPath((-0.3, 0), "A", size=1, prop=font_p),
            "V" : TextPath((-0.3, 0), "V", size=1, prop=font_p),
            "L" : TextPath((-0.3, 0), "L", size=1, prop=font_p),
            "I" : TextPath((-0.3, 0), "I", size=1, prop=font_p),
            "P" : TextPath((-0.3, 0), "P", size=1, prop=font_p),
            "W" : TextPath((-0.3, 0), "W", size=1, prop=font_p),
            "F" : TextPath((-0.3, 0), "F", size=1, prop=font_p),
            "M" : TextPath((-0.3, 0), "M", size=1, prop=font_p),
            "U" : TextPath((-0.3, 0), "U", size=1, prop=font_p),
            "O" : TextPath((-0.3, 0), "O", size=1, prop=font_p)}

aa_clr = {'C': 'limegreen', 'G': 'limegreen','S': 'limegreen','T': 'limegreen','Y': 'limegreen',
          'Q': 'purple','N': 'purple', 
          'K': 'blue','R': 'blue','H': 'blue', 
          'D': 'red','E': 'red',
          'A':'black','V':'black','L':'black','I':'black','P':'black','W':'black','F':'black','M':'black',
          'O':'yellow','U':'yellow'}


fig, ax = plt.subplots(figsize=(10,3))

x = 1
maxi = 0

for scores in aa_scores:
	y = 0
	for aa, score in scores:
		aa_info(aa, x,y, score, ax)
		y += score
	x += 1
	maxi = max(maxi, y)

plt.xticks(range(1,x))
plt.xlim((0, x)) 
plt.ylim((0, maxi)) 
plt.tight_layout()      
plt.show()
