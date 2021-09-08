#!/bin/bash

# Reference: 
# Garnier, J., Osguthorpe, D. J., & Robson, B. (1978).Analysis of the accuracy
# and implications of simple methods for predicting the secondary structure
# of globular proteins. Journal of molecular biology,120(1), 97–120.
# https://doi.org/10.1016/0022-2836(78)90297-8

if [[ -z "$1" ]]; then
	echo -e " Usage:\n\t$ bash GOR.sh <file.fasta>\t| Run with protein fasta record
\t$ bash GOR.sh --example\t\t| Run test on example data
\t$ bash GOR.sh --help\t\t| Show this help message and exit\n"


elif [[ $1 == "--help" ]]; then
	echo -e " Usage:\n\t$ bash GOR.sh <file.fasta>\t| Run with protein fasta record
\t$ bash GOR.sh --example\t\t| Run test on example data
\t$ bash GOR.sh --help\t\t| Show this help message and exit\n"

elif [[ $1 == "--example" ]]; then
	python bin/__main__.py example/example.fasta
else
	python bin/__main__.py $1

fi

