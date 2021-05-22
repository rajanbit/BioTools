# Chou-Fasman

Python implementation of Chou-Fasman algorithm for protein secondary structure prediction.

Chou-Fasman algorithm was the first protein structure prediction algorithm developed by Chou and Fasman in 1974. It pioneered a series of developments in protein structure prediction. This algorithm is based on a propensity table that was constructed by examining the x-ray-determined structures of 15 proteins containing 2473 amino acid residues and the number of occurrences of a given amino acid in the alpha-helix, beta-sheet, and coil was tabulated. The algorithm uses the propensity table and the set of empirical rules formulated by Chou and Fasman for protein secondary structure prediction.

>Requirements: Python 3.7.4

### ChouFas_predictor.py
Predict Alpha-Helical Segments, Beta-Sheet Segments and Beta-Turns within the protein sequence of a given protein fasta record
```
$ python ChouFas_predictor.py <fasta_file>
```

### helix_predictor.py
Predict helical segments within the protein sequence of a given protein fasta record
```
$ python helix_predictor.py <fasta_file>
```
### B-sheet_predictor.py
Predict beta-sheet segments within the protein sequence of a given protein fasta record
```
$ python B-sheet_predictor.py <fasta_file>
```
### turn_predictor.py
Predict beta-turn segments within the protein sequence of a given protein fasta record
```
$ python turn_predictor.py <fasta_file>
```

