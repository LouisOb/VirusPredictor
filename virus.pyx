from collections import Counter
from numpy import array
from Bio import SeqIO
import parmap

class virus:
	def __init__(self,virus_dna):
		if(type(virus_dna)==dict):
			self._DNA=virus_dna["RNA"]
			self._NAME=virus_dna["NAME"]
			self.ORFs=[]
			self.ORF_length=None
			self.na_length=None
		
	def ntuple(self,listing,N):
		if(self.na_length!=N):
			self.na_length=N
			Split_DNA=[self._DNA[i:i+N] for i in range(0,len(self._DNA),N)]
			Tuplecomb = dict(Counter(Split_DNA))
			for item in Tuplecomb:
				Tuplecomb[item]=Tuplecomb[item]/len(self._DNA)
			return Tuplecomb
	
	def ORF(self,length=12000):
		if(self.ORF_length!=length):
			self.ORF_length=length
			start=['A','T','G']
			stop=[['T','A','G'],['T','G','A'],['T','A','A']]
			starters=[[],[],[]]
			stoppers=[[],[],[]]
			for i in range(0,len(self._DNA)-2):
				if(self._DNA[i]==start[0] and self._DNA[i+1]==start[1] and self._DNA[i+2]==start[2]):
					starters[i%3].append(i)
				if(self._DNA[i]=='T'):
					for j in range(len(stop)):
						if(self._DNA[i+1]==stop[j][1] and self._DNA[i+2]==stop[j][2]):
							stoppers[i%3].append(i)
			self.ORFs=[]
			for k in range(len(starters)):
				stopdummy=0
				for i in range(len(starters[k])):
					if starters[k][i]<stopdummy:
						continue
					longest=[0]
					for stoppos in stoppers[k]:
						if(stoppos-starters[k][i]>2):
							longest=self._DNA[starters[k][i]:stoppos+3]
							stopdummy=stoppos
							if(len(longest)>=length):
								self.ORFs.append({"ORF":longest,"POS":starters[k][i],"ENDPOS":stoppos,"LEN":len(longest),"FRAME":k})	
							break

class virus_species:
	def __init__(self,virus_dna):
		if(type(virus_dna)==list):
			self.virus_dna_list=[]
			for element in virus_dna:
				self.virus_dna_list.append(virus(element))
	
	def ntuple(self,N):
		tuple_list = [self.virus_dna_list[i].ntuple(N) for i in range(len(self.virus_dna_list))]
		return tuple_list
	
	def na_count(self,na_comb):
		count=[]
		tuplelist=self.ntuple(len(na_comb[0]))
		for Tup in na_comb:
			for item in list(tuplelist):
				if Tup not in item.keys():
					tuplelist.remove(item)
					continue

		for item in tuplelist:
			for key in list(item.keys()):
				if key not in na_comb:
					del item[key]

		for Tup in na_comb:
			count.append(array([[tuplelist[i][Tup]] for i in range(len(tuplelist))]))
		return count


	def ORFs(self,length):
	
		ORF_list=[]
		processes=[]
		ORF_list = parmap.map(multi_ORF, self.virus_dna_list, length=length)
		return ORF_list

	def ORF_count(self,na_comb,length=12000):
		count=[]
		N=len(na_comb[0])
		ORFs=self.ORFs(length)
		if not len(ORFs):
			return array([])

		for listing in ORFs:
			for item in listing:
				Split_DNA=[item["ORF"][i:i+N] for i in range(0,len(item["ORF"]),N)]
				Tuplecomb = dict(Counter(Split_DNA))
				for item1 in Tuplecomb:
					Tuplecomb[item1]=Tuplecomb[item1]/len(item["ORF"])
				count.append(Tuplecomb)
		for Tup in na_comb:
			for item in count:
				if Tup not in item.keys():
					item.update({Tup:0})

		for item in count:
			for key in list(item.keys()):
				if key not in na_comb:
					del item[key]
		result=[]
		for Tup in na_comb:
			result.append([count[i][Tup] for i in range(len(count))])

		return array(result)

def readfile(filename):
	DNA=[]
	with open(filename,"r") as fastafile:
		for sequence in SeqIO.parse(fastafile,"fasta"):
			DNA.append({"RNA":sequence.seq,"NAME":sequence.name})			
	return DNA;

def multi_ORF(item,length):
			item.ORF(length)
			return item.ORFs
