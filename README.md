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

## Clustering.py
Contains the funtions to extract features from nucleobase combinations of virus genome sequences. The combination of nucleobased used as features is passed to initialize the clustering object.
```python
import clustering

#Introduce Tuple that are searched for in the sequence
	NA = ["A","C","G","T"]
	Tup = []
	for na1 in NA:
		for na2 in NA:
			for na3 in NA:
						Tup.append(na1+na2+na3)

	viruslist=clustering(Tup)
	viruslist.add(["SARS_COV_2_new_new.fasta","SARS-CoV-1.fasta.txt","MERS.fasta","Bovine_coronavirus.fasta","Camel_alphacoronavirus.fasta","Duck_coronavirus.fasta"],["SARS CoV 2","SARS","MERS","Bovine corona virus","Camel alpha corona virus","Duck corona virus"])
	viruslist.ORF_tuples(11000)
	viruslist.PCA(2,True)
	viruslist.plot("Coronafamily_test.png")
```
