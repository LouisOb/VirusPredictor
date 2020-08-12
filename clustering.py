#!/usr/bin/python3
import sklearn
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pickle as pk
import numpy as np
from virus import virus, virus_species, readfile

pca = sklearn.decomposition.PCA(n_components=2) 

#load DNA of viruses
DNA_MERS=readfile("MERS.fasta")
DNA_SARS_COV_1=readfile("SARS-CoV-1.fasta.txt")
DNA_SARS_COV_2=readfile("SARS_COV_2.fasta")
DNA_HEPA=readfile("Hepatitis_C_1b.fasta")
DNA_HERPES=readfile("HERPES.fasta")

#Introduce Tuple that are searched for in the sequence
NA = ["A","C","G","T"]
Tup = []
for na1 in NA:
	for na2 in NA:
		for na3 in NA:
					Tup.append(na1+na2+na3)

#Setup the virus DNA Objects

HERPES_LIST=virus_species(DNA_HERPES)
HERPES_LIST_plot=HERPES_LIST.na_count(Tup)


#HERPES_comb = np.concatenate(HERPES_LIST_plot[:],axis=1)


MERS_list=virus_species(DNA_MERS)
MERS_list_plot = MERS_list.na_count(Tup)

#MERS_comb = np.concatenate(MERS_list_plot[:],axis=1)

SARS_1_list=virus_species(DNA_SARS_COV_1)
SARS_list_plot = SARS_1_list.na_count(Tup)

#SARS_1_comb = np.concatenate(SARS_list_plot[:],axis=1)

SARS_2_list=virus_species(DNA_SARS_COV_2)
SARS_2_list_plot = SARS_2_list.na_count(Tup)

#SARS_2_comb = np.concatenate(SARS_2_list_plot[:],axis=1)

HEPA_list = virus_species(DNA_HEPA)
HEPA_list_plot=HEPA_list.na_count(Tup)

#HEPA_comb = np.concatenate(HEPA_list_plot[:],axis=1)

#Virus sequence DNA PCA with selected Tuples
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

#Find ORFs in Virus DNA sequence and count Tuple occurence

SARS_2_ORFs=SARS_2_list.ORF_count(Tup)
SARS_1_ORFs=SARS_1_list.ORF_count(Tup)
MERS_ORFs=MERS_list.ORF_count(Tup)
#HERPES_ORFs=HERPES_LIST.ORF_count(Tup)
#HEPA_ORFs=HEPA_list.ORF_count(Tup)

colour_mers=['blue' for i in range(len(MERS_ORFs.T))]
colour_sars_1=['red' for i in range(len(SARS_1_ORFs.T))]
colour_sars_2=['green' for i in range(len(SARS_2_ORFs.T))]
#colour_hepa=['purple' for i in range(len(HEPA_ORFs.T))]
#colour_herpes=['black' for i in range(len(HERPES_ORFs.T))]

colour=colour_mers+colour_sars_1+colour_sars_2#+colour_herpes+colour_hepa

virus_comb = np.concatenate((MERS_ORFs,SARS_1_ORFs,SARS_2_ORFs),axis=1)#,HERPES_ORFs,HEPA_ORFs),axis=1)


#PCA with ORF Tuple count data

#Transformation_Mat = pk.load(open("3tuple_pca.pk","rb"))
#DNA_ORFs=np.matmul(Transformation_Mat,virus_comb)


virus_comb = pca.fit_transform(virus_comb.T)
pk.dump(pca.components_,open("3tuple_pca_new.pk","wb"))


#Histograms with the Tuple distribution
'''
for i in range(len(Tup)):
		plt.hist([SARS_2_ORFs[i][j] for j in range(len(SARS_2_ORFs[i]))],bins=20,alpha=.5,facecolor="green",label="SARS-CoV 2")
		plt.hist([SARS_1_ORFs[i][j] for j in range(len(SARS_1_ORFs[i]))],bins=20,alpha=.5,facecolor="red",label="SARS")
		plt.hist([MERS_ORFs[i][j] for j in range(len(MERS_ORFs[i]))],bins=20,alpha=.5,facecolor="blue",label="MERS")
#		plt.hist(virus_comb[i],10,color=["blue","red","green"])
		plt.xlabel(Tup[i])
		plt.legend()
		plt.savefig("pics/"+Tup[i]+".png",pad_inches=0.1,bbox_inches='tight')
#		plt.show()
		plt.close()
'''

plt.scatter(virus_comb.T[0],virus_comb.T[1],c=colour)
plt.title(str(Tup)+" PCA Feature Space")
plt.show()








