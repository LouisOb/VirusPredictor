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
	"""
	The class with the means to extract features based on nucleobase combinations and PCA.

	A clustering class object is able to load fasta files and extract features based on a
	preselected set of nucleobase combinations. The features are extracted by counting the
	occurences of a nucleobase combination and dividing it by the size of the whole
	sequence.
	PCA analysis can then be applied to the vector containing the normalized count to
	reduce the dimensionality to a selected value.
	Attributes:	
		DNA_set: 	A list containing lists for each virus file supplied containing
				 	dictionaries with the DNA sequences of the viruses.
		viruslist:	A list containing objects of the class virus.
		colors:		A list containing plotting colors.
		plot_colors:The colors used when calling the plotting function.
		tuples:		The list containing the nucleobase combinations which are counted.
		tuple_count:The list containing the count of the nucleobase combinations for the
					supplied virus sequences.
		orf_count:	The list containing the count of nucleobase combinations from ORFs of
					virus sequences.
		plot_list:	The list including the two list to be plotted.
		names:		The list containing the names of the fasta files or a supplied name.
		mode:		A string which saves the operation used on the set. Possible values are
					"orf" if the extraction has been executed of an ORF or "tuples" if the whole
					sequence has been analysed.
	
	"""
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
		"""
		The Function which adds virus files to the set.
		
		This function loads fasta formated files and adds them to the list of virus sequences.
		Args:
			filename: Either a string or a list containing the name or names of the fasta
					  file/s with the virus sequence/s.
			name:	  A string or list with the name/s to be assigned to the virus 
					  sequence/s from the fasta files/s.
		Returns:
			Nothing
		"""
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
		"""
		The Function which counts the nucleobase combinations.
		
		This function counts the nucleobase combinations in the whole sequence of all
		virus sequences contained in the object viruslist. The extracted features are
		stored in the list tuple_count.
		Args:
			None
		Returns:
			Nothing
		"""
		
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
		"""
		The Function which counts the nucleobase combinations in an ORF.
		
		This function counts the nucleobase combinations in ORFs of a specified size
		of all virus sequences contained in the object viruslist. The extracted 
		features are stored in the list tuple_count.
		Args:
			lenght: The minimum length of the ORF sequence used to extract features.
		Returns:
			Nothing
		"""
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
		"""
		The Function which carries out PCA.
		
		This function uses PCA on the list tuple_count and saves the results to the
		plot_list.
		Args:
			components: The selected dimensionality the result has to have after PCA.
			use_old:	Use a PCA Matrix saved in an older run to transform the set.
		Returns:
			Nothing
		"""
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
		"""
		The Function which plots two dimensional data.
		
		This function uses the matplotlib to plot the two dimensional data contained
		in the list plot_list.
		Args:
			filename: The filename which the image will be saved to.
			title:	  The title of the matplotlib plot.
		Returns:
			Nothing
		"""
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
		"""
		The Function which loads pickle files.
	
		The function used to load pickle files containing information or
		results saved in previous runs.

		Args:
			filename: The name of the pickle file.
		Returns:
			Nothing
		"""
		return pk.load(open("save/"+filename,"rb"))
		
	def save(self):
		"""
		The Function which saves to pickle files.
	
		The function used to save extraced results to pickle files.

		Args:
			None
		Returns:
			Nothing
		"""
		for index in range(len(self.names)):
			if(self.tuple_count==[]):
				dump_array = [[],self.orf_count[index]]
			elif(self.orf_count==[]):
				dump_array = [self.tuple_count[index],[]]
			else:
				dump_array = [self.tuple_count[index],self.orf_count[index]]
			pk.dump(dump_array,open("save/"+self.names[index][1]+".pk","wb"))
	
	def clean(self):
		"""
		A Function which clears some lists
	
		This function clears some lists which enabels a new feature extraction.

		Args:
			None
		Returns:
			Nothing
		"""
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
	testlist.add(["SARS_COV_2_new_new.fasta","SARS-CoV-1.fasta.txt","MERS.fasta","Bovine_coronavirus.fasta","Camel_alphacoronavirus.fasta","Duck_coronavirus.fasta"],["SARS CoV 2","SARS","MERS","Bovine corona virus","Camel alpha corona virus","Duck corona virus"])
	testlist.ORF_tuples(11000)
	testlist.PCA(2,True)
	testlist.plot("Coronafamily_test.png")









