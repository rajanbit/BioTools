import sys

inpt = open(sys.argv[1],"r")
inpt_rl = inpt.readlines()
outpt = open(sys.argv[1]+".tsv","w+")
sample = ""
head = "Sample\t"
row = ""

for line in inpt_rl:
	data = line.strip().split("\t")
	sample = (data[0])+"\t"
	head += data [1]+"\t"
	row += data[2]+"\t"

outpt.write(head+"\n")
outpt.write(sample+row+"\n")
