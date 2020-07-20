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
Tup=['GTTG','CAAG','AATT','ACAA','CAAC','TTTA','GTTA','TCAT','TTAA']


HERPES_LIST=virus_species(DNA_HERPES)
HERPES_LIST_plot=HERPES_LIST.na_count(Tup)

HERPES_comb = np.concatenate(HERPES_LIST_plot[:],axis=1)


MERS_list=virus_species(DNA_MERS)
MERS_list_plot = MERS_list.na_count(Tup)

MERS_comb = np.concatenate(MERS_list_plot[:],axis=1)

SARS_1_list=virus_species(DNA_SARS_COV_1)
SARS_list_plot = SARS_1_list.na_count(Tup)

SARS_1_comb = np.concatenate(SARS_list_plot[:],axis=1)

SARS_2_list=virus_species(DNA_SARS_COV_2)
SARS_2_list_plot = SARS_2_list.na_count(Tup)

SARS_2_comb = np.concatenate(SARS_2_list_plot[:],axis=1)

HEPA_list = virus_species(DNA_HEPA)
HEPA_list_plot=HEPA_list.na_count(Tup)

HEPA_comb = np.concatenate(HEPA_list_plot[:],axis=1)
'''
colour_mers=['blue' for i in range(len(MERS_comb))]
colour_sars_1=['red' for i in range(len(SARS_1_comb))]
colour_sars_2=['green' for i in range(len(SARS_2_comb))]
colour_hepa=['purple' for i in range(len(HEPA_comb))]
colour_herpes=['black' for i in range(len(HERPES_comb))]
colour=colour_mers+colour_sars_1+colour_sars_2+colour_hepa#+colour_herpes


virus_comb = np.concatenate((MERS_comb,SARS_1_comb,SARS_2_comb,HEPA_comb),axis=0)
virus_comb=pca.fit_transform(virus_comb)'''
#plt.scatter(virus_comb.T[0],virus_comb.T[1],c=colour)
#plt.legend()
#plt.show()

SARS_2_ORFs=SARS_2_list.ORF_count(Tup)
SARS_1_ORFs=SARS_1_list.ORF_count(Tup)
MERS_ORFs=MERS_list.ORF_count(Tup)
#HERPES_ORFs=HERPES_LIST.ORF_count(Tup)
#HEPA_ORFs=HEPA_list.ORF_count(Tup)

#pk.dump(MERS_ORFs,open("MERS_ORFs.pk","wb"))
#pk.dump(SARS_1_ORFs,open("SARS_ORFs.pk","wb"))
#pk.dump(SARS_2_ORFs,open("SARS_2_ORFs.pk","wb"))
#SARS_2_ORFs=pk.load(open("SARS_2_ORFs.pk","rb"))
#SARS_1_ORFs=pk.load(open("SARS_ORFs.pk","rb"))
#MERS_ORFs=pk.load(open("MERS_ORFs.pk","rb"))

SARS_1_comb = np.concatenate(SARS_1_ORFs[:],axis=1)
SARS_2_comb = np.concatenate(SARS_2_ORFs[:],axis=1)
MERS_comb = np.concatenate(MERS_ORFs[:],axis=1)
#HERPES_comb = np.concatenate(HERPES_ORFs[:],axis=1)
#HEPA_comb = np.concatenate(HEPA_ORFs[:],axis=1)

colour_mers=['blue' for i in range(len(MERS_comb))]
colour_sars_1=['red' for i in range(len(SARS_1_comb))]
colour_sars_2=['green' for i in range(len(SARS_2_comb))]
colour_hepa=['purple' for i in range(len(HEPA_comb))]
colour_herpes=['black' for i in range(len(HERPES_comb))]


colour=colour_mers+colour_sars_1+colour_sars_2

virus_comb = np.concatenate((MERS_comb,SARS_1_comb,SARS_2_comb),axis=0)
virus_comb = pca.fit_transform(virus_comb)
#plt.scatter(SARS_2_ORFs[0],SARS_2_ORFs[1],c='green',marker="x")
#plt.scatter(SARS_1_ORFs[0],SARS_1_ORFs[1],c='red',marker="+")
#plt.scatter(MERS_ORFs[0],MERS_ORFs[1],c='blue',marker="^")


'''
for i in range(len(Tup)):
		fig=plt.figure()
		ax=fig.add_axes([0,0,1,1])
		ax.hist(SARS_2_ORFs[i],label="SARS-CoV 2")
		ax.hist(SARS_1_ORFs[i],label="SARS")
		ax.hist(MERS_ORFs[i],label="MERS")
		ax.set_xlabel(Tup[i])
#		ax.set_ylabel(Tup[j])
		ax.legend()
		plt.savefig("pics/"+Tup[i]+".png",pad_inches=0.1,bbox_inches='tight')
		plt.close(fig)
'''

plt.scatter(virus_comb.T[0],virus_comb.T[1],c=colour)
plt.title(str(Tup)+" PCA Feature Space")
plt.show()






