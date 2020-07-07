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
	first=True
	for line in fileo:
		if(first==False and line[0]=='>'):
			dummy.update({"RNA":dummystring.replace('\n',"")})
			DNA.append(dummy)
			dummy = {}
			dummy.update({"NAME":line[1:].replace('\n',"")})
			dummystring=""
			continue;
		if(line[0]==">" and first):
			dummy.update({"NAME":line[1:].replace('\n',"")})
			first=False
			continue;
		dummystring=dummystring+line
	return DNA;


pca = sklearn.decomposition.PCA(n_components=2) 

DNA_MERS=readfile("MERS.fasta")
DNA_SARS_COV_1=readfile("SARS-CoV-1.fasta.txt")
DNA_SARS_COV_2=readfile("SARS_COV_2.fasta")
DNA_HEPA=readfile("Hepatitis_C_1b.fasta")
DNA_HERPES=readfile("HERPES.fasta")

NA = ["A","C","G","T"]
Tup = []
for na1 in NA:
	for na2 in NA:
		for na3 in NA:
			for na4 in NA:
				Tup.append(na1+na2+na3+na4)	


HERPES_LIST=virus_species(DNA_HERPES)
HERPES_LIST=HERPES_LIST.na_count(Tup)

HERPES_comb = np.concatenate(HERPES_LIST[:],axis=1)


MERS_list=virus_species(DNA_MERS)
MERS_list = MERS_list.na_count(Tup)

MERS_comb = np.concatenate(MERS_list[:],axis=1)

SARS_1_list=virus_species(DNA_SARS_COV_1)
SARS_list = SARS_1_list.na_count(Tup)

SARS_1_comb = np.concatenate(SARS_list[:],axis=1)

SARS_2_list=virus_species(DNA_SARS_COV_2)
SARS_2_list = SARS_2_list.na_count(Tup)

SARS_2_comb = np.concatenate(SARS_2_list[:],axis=1)

HEPA_list = virus_species(DNA_HEPA)
HEPA_list=HEPA_list.na_count(Tup)

HEPA_comb = np.concatenate(HEPA_list[:],axis=1)

colour_mers=['blue' for i in range(len(MERS_comb))]
colour_sars_1=['red' for i in range(len(SARS_1_comb))]
colour_sars_2=['green' for i in range(len(SARS_2_comb))]
colour_hepa=['purple' for i in range(len(HEPA_comb))]
colour_herpes=['black' for i in range(len(HERPES_comb))]
colour=colour_mers+colour_sars_1+colour_sars_2+colour_hepa+colour_herpes


virus_comb = np.concatenate((MERS_comb,SARS_1_comb,SARS_2_comb,HEPA_comb,HERPES_comb),axis=0)
virus_comb=pca.fit_transform(virus_comb)

plt.scatter(virus_comb[:,0],virus_comb[:,1],c=colour,label="PCA of Featurespace")
plt.legend()
plt.show()








