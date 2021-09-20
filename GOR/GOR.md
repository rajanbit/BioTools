# GOR

GOR (Garnier–Osguthorpe–Robson) is a protein secondary structure prediction method developed in late 1970s shortly after Chou–Fasman method.
This method is also based on the probability tables constructed examining the information on protein structures derieved from X-ray crystallography
but it also takes into consideration the probability values of immediate neighbors. GOR method predict alpha helix, beta sheet, turn, or random coil
at each position based on 17-amino-acid sequence windows.

### Requirements:

> OS Linux
>
> Python 3.7.4

### Usage:

```

 Usage:
	$ bash GOR.sh <file.fasta>	| Run with protein fasta record
	$ bash GOR.sh --example		| Run test on example data
	$ bash GOR.sh --help		| Show this help message and exit

```
