''' \
Usage:
   python b.py -n <file.fasta> -o <output.tsv> -c <coverage>'''

import sys

...

inputs = sys.argv

if '-n' not in inputs and '-o' not in inputs:
	print (__doc__)

else:
	f_in = inputs[inputs.index('-n') + 1]
	f_out = inputs[inputs.index('-o') + 1]
	cov = int(inputs[inputs.index("-c") + 1])

	def gc_analyzer():

		f_open = open(f_in, "r")
		f_read = f_open.readlines()
		f_write = open(f_out, "w+")

		head = ""
		seq = ""
		for line in f_read:
			if ">" in line[0:2]:
				head = line
			elif ">" not in line[0:2]:
				seq += line.strip()
		for i in range(0, len(seq), cov):
			seq_frag = seq[i:i+cov]
			gc = ("%.2f" %(((seq_frag.count("G")+seq_frag.count("C"))/len(seq_frag))*100))
			f_write.write(str((i+cov)/1000000)+"\t"+str(gc)+"\n")
	gc_analyzer()

...
