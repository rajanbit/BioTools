import sys

fna_fl = sys.argv[1]
gff_fl = sys.argv[2]

fna_op = open(fna_fl, "r")
gff_op = open(gff_fl, "r")

fna = fna_op.readlines()
seq = ""
for line in fna:
	if line[0:1] != ">":
		seq += line.strip()

gff = gff_op.readlines()
for line in gff:
	if line[0:1] != "#":
		data = line.split(";")
		if "gbkey=CDS" in data: 
			head = ">"+data[0].split("\t")[0]+"_"+data[4][6:]+"_"+data[3][5:]
			if "partial" not in head and "incomplete" not in head and "frameshift" not in head and len(head) <= 50 and "stop" not in head:
				head = head+" ["+data[-4]+"]"+" ["+data[2].split(",")[1]+"]"+" ["+data[-3]+"]"+" ["+data[-2]+"]"
				if data[0].split("\t")[6] == "+":
					cds_seq = seq[int(data[0].split("\t")[3])-1:int(data[0].split("\t")[4])]
					temp_fasta = ""
					for i in range(0, len(cds_seq), 70):
						temp_fasta += cds_seq[i:i+70]+"\n"
					print(head)
					print(temp_fasta)
				
