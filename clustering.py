#!/usr/bin/python3
import sklearn
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pickle as pk
import numpy as np
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


pca = sklearn.decomposition.PCA(n_components=2) 

DNA_MERS=readfile("MERS.fasta.txt")
DNA_SARS_COV_1=readfile("SARS-CoV-1.fasta.txt")
DNA_SARS_COV_2=readfile("SARS-CoV-2.fasta.txt")

Tup = ["GCA","GCC","AAA","TCC","ACT","GTT","GAC","GTA"]

MERS_list=virus_species(DNA_MERS)
MERS_1,MERS_2,MERS_3,MERS_4,MERS_5,MERS_6,MERS_7,MERS_8 = MERS_list.na_count(Tup)

MERS_comb = np.concatenate((MERS_1,MERS_2,MERS_3,MERS_4,MERS_5,MERS_6,MERS_7,MERS_8),axis=1)

SARS_1_list=virus_species(DNA_SARS_COV_1)
SARS_1_1,SARS_1_2,SARS_1_3,SARS_1_4,SARS_1_5,SARS_1_6,SARS_1_7,SARS_1_8 = SARS_1_list.na_count(Tup)

SARS_1_comb = np.concatenate((SARS_1_1,SARS_1_2,SARS_1_3,SARS_1_4,SARS_1_5,SARS_1_6,SARS_1_7,SARS_1_8),axis=1)

SARS_2_list=virus_species(DNA_SARS_COV_2)
SARS_2_1,SARS_2_2,SARS_2_3,SARS_2_4,SARS_2_5,SARS_2_6,SARS_2_7,SARS_2_8 = SARS_2_list.na_count(Tup)

SARS_2_comb = np.concatenate((SARS_2_1,SARS_2_2,SARS_2_3,SARS_2_4,SARS_2_5,SARS_2_6,SARS_2_7,SARS_2_8),axis=1)

colour_mers=['blue' for i in range(len(MERS_comb))]
colour_sars_1=['red' for i in range(len(SARS_1_comb))]
colour_sars_2=['green' for i in range(len(SARS_2_comb))]
colour=colour_mers+colour_sars_1+colour_sars_2


virus_comb = np.concatenate((MERS_comb,SARS_1_comb,SARS_2_comb),axis=0)
virus_comb=pca.fit_transform(virus_comb)

plt.scatter(virus_comb[:,0],virus_comb[:,1],c=colour,label="PCA of Featurespace")
plt.legend()
plt.show()








