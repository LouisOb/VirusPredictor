# Virus Classification

A project which engages with the topic of virus classification.
The classification is based on RNA sequences as a result of a sequencing tool.

## Virus.pyx
Contains the main functions used to count Nucleobase tuples or search for ORFs in sequences.
The virus and virus_species classes are used to achieve this.
Read in sequences are passed for initialization of the virus_species object, which itself will initialize virus objects.
The virus species object can then be used to return Tuple counts or ORFs.

```python
from virus import virus, virus_species, readfile

DNA_SARS_COV_2=readfile("SARS_COV_2.fasta")
#Introduce Tuple that are searched for in the sequence
NA = ["A","C","G","T"]
Tup = []
for na1 in NA:
	for na2 in NA:
		for na3 in NA:
					Tup.append(na1+na2+na3)

SARS_2_list=virus_species(DNA_SARS_COV_2)
SARS_2_list_plot = SARS_2_list.na_count(Tup)
```

## Fingerprint.py
Contains the functions to identify a common pattern of nucleobases and aminoacids after translation in a virus DNA sequence.
A fasta file which contains the Virus DNA sequences has to be passed to the functions.

```python
import Fingerprint

Fingerprint.ORF_scan("SARS-CoV-1.fasta.txt")
Fingerprint.test_identity("SARS-CoV-1.fasta.txt")
```

It will then print the percentage of the sequences in the file containing the discovered pattern.
