#!/usr/bin/python3
import sklearn
import matplotlib.pyplot as plt
import pickle as pk
from virus import virus, virus_species

def readfile(filename):
	fileo = open(filename)
	DNA=[]
	dummy = {}
	dummystring=""
	for line in fileo:
		if(line[0]=='\n'):
			dummy.update({"RNA":dummystring.replace('\n',"")})
			DNA.append(dummy)
			dummy = {}
			dummystring=""
			continue;
		if(line[0]==">"):
			dummy.update({"NAME":line[1:].replace('\n',"")})
			continue;
		dummystring=dummystring+line
	return DNA;



DNA_MERS=readfile("MERS.fasta.txt")
DNA_SARS_COV_1=readfile("SARS-CoV-1.fasta.txt")
DNA_SARS_COV_2=readfile("SARS-CoV-2.fasta.txt")

MERS_list=virus_species(DNA_MERS)
MERS_x,MERS_y = MERS_list.na_count("AG","GT")

plt.scatter(MERS_x,MERS_y,c='blue',label='MERS')
#plt.scatter(SARS_1_x,SARS_1_y,c='red',label='SARS_COV_1')
#plt.scatter(SARS_2_x,SARS_2_y,c='green',label='SARS_COV_2')
plt.legend()
plt.xlabel("AG")
plt.ylabel("GT")
plt.show()








