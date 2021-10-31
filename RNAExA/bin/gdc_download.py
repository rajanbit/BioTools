#!/usr/bin/python
#!python

import sys
import os

file1 = sys.argv[1]
gdc_manifest_file = open(file1, "r")
gdc_data = gdc_manifest_file.readlines()
log = open("temp/logfile.txt", "w")
for data in gdc_data[1:]:
	data_ls = data.split("\t")
	gdc_id = data_ls[0]
	os.system("wget -q https://api.gdc.cancer.gov/data/"+gdc_id+" -O HTSeq_Counts/"+gdc_id+".htseq.counts.gz")
	log = open("temp/logfile.txt", "a")
	log.write(gdc_id+"\n")
log.close()
gdc_manifest_file.close()

# python gdc_download.py <gdc_manifest.txt>
