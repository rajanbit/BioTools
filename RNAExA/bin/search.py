#!/usr/bin/python
#!python

import sys
import os

imput_tab = open(sys.argv[1], "r")
htseq_count = open(sys.argv[2], "r")
out_fname = "temp/"+sys.argv[2]
out_fname = out_fname.split("/")[-1]
out_file = open("temp/"+out_fname+".out", "w")
imput_tab_rl = imput_tab.readlines()
htseq_count_rl = htseq_count.readlines()
for line in htseq_count_rl:
	ensemble_id = line.split("\t")[0]
	for data in imput_tab_rl:
		if ensemble_id.split(".",1)[0] == data.split("\t")[0].split(".",1)[0]:
			out_file.write(line.strip()+"\t"+data.split("\t")[1]+"\t"+data.split("\t")[2]+"\n")


imput_tab.close()
htseq_count.close()
out_file.close()
