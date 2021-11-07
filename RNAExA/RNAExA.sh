#!/bin/bash

if [[ -z "$1" ]]; then
	echo -e ""


elif [[ $1 == "--help" ]]; then
	echo -e ""

elif [[ $1 == "--download"|| $1 == "-d" ]]; then
	python bin/gdc_download.py Manifest/*
	gunzip HTSeq_Counts/*.gz

elif [[ $1 == "--example"|| $1 == "-e" ]];then
	cp example_data/gdc_manifest.txt Manifest/
	cp example_data/input.tab input.tab
	python bin/gdc_download.py Manifest/*
	gunzip HTSeq_Counts/*.gz
	for i in HTSeq_Counts/*.htseq.counts; do python bin/search.py input.tab $i;done
	for i in temp/*.out; do awk '{print FILENAME"\t"$3"\t"$2}' $i |sed "s/.htseq.counts.out//g" | sed "s/[temp/]//g"> $i.a;done
	for i in temp/*.out.a; do python bin/make_tsv.py $i;done
	cat temp/*.tsv > temp_out.tsv
	awk '! a[$0]++' temp_out.tsv > output.tsv
	rm temp_out.tsv
	Rscript bin/boxplot.R output.tsv

elif [[ $1 == "--search"|| $1 == "-s" ]];then
	rm temp/*
	for i in HTSeq_Counts/*.htseq.counts; do python bin/search.py input.tab $i;done
	for i in temp/*.out; do awk '{print FILENAME"\t"$3"\t"$2}' $i |sed "s/.htseq.counts.out//g" | sed "s/[temp/]//g"> $i.a;done
	for i in temp/*.out.a; do python bin/make_tsv.py $i;done
	cat temp/*.tsv > temp_out.tsv
	awk '! a[$0]++' temp_out.tsv > output.tsv
	rm temp_out.tsv

elif [[ $1 == "--plot"|| $1 == "-p" ]];then
	Rscript bin/boxplot.R output.tsv

elif [[ $1 == "--clear"|| $1 == "-c" ]];then
	rm temp/*
	rm HTSeq_Counts/*
	rm Manifest/*
	rm output.tsv
	rm input.tab
	rm boxplot.png

fi
