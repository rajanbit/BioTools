# GC_Classifier
Check the genomic GC profile based similarity between the two fasta record (Reference genome and Test genome) and plot the data to show the relationship between the two genome

### Requirements:
> OS Linux

> Python 3.7.4

> R version 3.5.1 (2018-07-02) -- "Feather Spray"

> Library ggplot2

### Usage:
Paste the FASTA files (file1.fasta & file2.fasta) within the bin folder and run the script classifier.sh
```bash
$ bash classifier.sh <file1.fasta> <file2.fasta>
```
### Output:
The program generate a result.png file within the bin folder, in which the GC profile of both the samples are plotted
