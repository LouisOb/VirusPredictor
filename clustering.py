#!/usr/bin/python3
import pyximport; pyximport.install()
import sklearn
from sklearn.decomposition import PCA
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import pickle as pk
import numpy as np
from virus import virus, virus_species, readfile
from os import listdir

class clustering:
	DNA_set		= []
	viruslist	= None
	colors		= ["green","red","blue","black","purple","orange","pink","yellow","darkgoldenrod","deepskyblue","grey","lime","cyan"]
	plot_colors	= []
	tuples		= None
	tuple_count	= []
	orf_count	= []
	plot_list	= []
	names 		= []
	mode		= None

	def __init__(self,tuples_list):
		self.tuples=tuples_list

	def add(self,filename,name=[]):
		if(type(filename)==str):
			if(name==[]):
				self.names.append([filename,filename])
			else:
				self.names.append([name,filename])
			DNA	= readfile("resources/"+filename)
			self.DNA_set.append(DNA)
		if(type(filename)==list):
			if(name==[]):
				for entry in filename:
					self.names.append([entry,entry])
			else:
				for index in range(len(name)):
					self.names.append([name[index],filename[index]])
			for entry in filename:
				DNA	= readfile("resources/"+entry)
				self.DNA_set.append(DNA)

	def tuples(self):
		self.clean()
		self.mode = "tuples"
		for listing in self.DNA_set:
			self.viruslist.append(virus_species(listing))
		for index in range(len(self.viruslist)):
			if(self.names[index][1]+".pk" in listdir("./save")):
				check = self.load(self.names[index][1]+".pk")[0]
				if( check.size > 0): 	
					self.tuple_count.append(check)
					continue
			self.tuple_count.append(item.na_count(self.tuples))
		self.save()
		color_index=0
		for entry in self.tuple_count:
			for i in range(len(entry)):
				self.plot_colors.append(self.colors[color_index])
			color_index+=1
	
	def ORF_tuples(self,length=12000):
		self.clean()
		self.mode = "orf"
		for listing in self.DNA_set:
			self.viruslist.append(virus_species(listing))
		for index in range(len(self.viruslist)):
			if(self.names[index][1]+".pk" in listdir("./save")):
				check = self.load(self.names[index][1]+".pk")[1]
				if( check.size > 0): 
					self.orf_count.append(check)
					continue
			self.orf_count.append(self.viruslist[index].ORF_count(self.tuples,length))
		self.save()
		color_index=0
		for entry in self.orf_count:
			for i in range(len(entry[0])):
				self.plot_colors.append(self.colors[color_index])
			color_index+=1

	def PCA(self,components,use_old=False):
		if(self.mode=="tuples"):
			pca_list=np.hstack(self.tuple_count)
		else:
			pca_list=np.hstack(self.orf_count)
		if(use_old):
			Transformation_Mat 	= pk.load(open("Mat_pca.pk","rb"))
			self.plot_list		= np.matmul(Transformation_Mat,pca_list).T
		else:
			pca 			= sklearn.decomposition.PCA(n_components=components)
			self.plot_list 	= pca.fit_transform(pca_list.T)
			pk.dump(pca.components_,open("Mat_pca.pk","wb"))

	def plot(self,filename,title=""):
		if(not len(self.plot_list)):
			print("[Warning] No plot list availiable for object")
			return
		legend_elements = []
		for index in range(len(self.names)):
			legend_elements.append(Line2D([0], [0], marker='.', color="w", label=self.names[index][0], markerfacecolor=self.colors[index], markersize=11))
		plt.scatter(self.plot_list.T[0],self.plot_list.T[1],c=self.plot_colors,marker='.')
		plt.legend(handles=legend_elements, loc='best', borderaxespad=0.)
		plt.title(title)

		plt.gcf().set_size_inches(10, plt.gcf().get_size_inches()[1])
		plt.xlabel("PCA feature combination 0")
		plt.ylabel("PCA feature combination 1")
		plt.savefig("pics/"+filename,bbox_inches='tight',dpi=250)
		plt.close()

	def load(self,filename):
		return pk.load(open("save/"+filename,"rb"))
		
	def save(self):
		for index in range(len(self.names)):
			if(self.tuple_count==[]):
				dump_array = [[],self.orf_count[index]]
			elif(self.orf_count==[]):
				dump_array = [self.tuple_count[index],[]]
			else:
				dump_array = [self.tuple_count[index],self.orf_count[index]]
			pk.dump(dump_array,open("save/"+self.names[index][1]+".pk","wb"))
	
	def clean(self):
		self.viruslist		= []
		self.plot_colors	= []
		self.tuple_count	= []


if __name__ == "__main__":

	#Introduce Tuple that are searched for in the sequence
	NA = ["A","C","G","T"]
	Tup = []
	for na1 in NA:
		for na2 in NA:
			for na3 in NA:
						Tup.append(na1+na2+na3)

	testlist=clustering(Tup)
#	testlist.add("Duck_coronavirus.fasta")
	testlist.add(["SARS_COV_2_new_new.fasta","SARS-CoV-1.fasta.txt","MERS.fasta","Bovine_coronavirus.fasta","Camel_alphacoronavirus.fasta","Duck_coronavirus.fasta"],["SARS CoV 2","SARS","MERS","Bovine corona virus","Camel alpha corona virus","Duck corona virus"])
#	testlist.add(["SARS_COV_2_new.fasta","SARS_COV_1.fasta","MERS.fasta"],["SARS CoV 2","SARS","MERS"])
#	testlist.add(["Alpha_Herpes_virus.fasta","Beta_Herpes_virus.fasta","Gamma_Herpes_virus.fasta"],["Alpha herpes virus","Beta herpes virus","Gamma herpes virus"])
#	testlist.add(["Mastadeno_A.fasta","Mastadeno_B.fasta","Mastadeno_C.fasta","Mastadeno_D.fasta","Mastadeno_E.fasta","Mastadeno_F.fasta"],["Mastadeno A virus", "Mastadeno B virus", "Mastadeno C virus", "Mastadeno D virus", "Mastadeno E virus", "Mastadeno F virus"])
#	testlist.add(["Vaccinia_virus.fasta","Variola_virus.fasta","Cowpox_virus.fasta"],["Vaccinia virus", "Variola virus", "Cowpox virus"])
	testlist.ORF_tuples(11000)
	testlist.PCA(2,True)
	testlist.plot("SARSfamily_ORF11000.png")
#	testlist.plot("herpes_corona_family.png")







