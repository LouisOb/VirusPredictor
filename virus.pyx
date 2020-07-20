from collections import Counter
from numpy import array

class virus:
	def __init__(self,virus_dna):
		if(type(virus_dna)==dict):
			self._DNA=virus_dna["RNA"]
			self._NAME=virus_dna["NAME"]
			self.ORFs=None
		
	def ntuple(self,listing,N):
		Split_DNA=[listing[i:i+N] for i in range(0,len(listing),N)]
		Tuplecomb = dict(Counter(Split_DNA))
		for item in Tuplecomb:
			Tuplecomb[item]=Tuplecomb[item]/len(listing)
		return Tuplecomb
	
	def ORF(self,length=10000):
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
		tuple_list = [self.virus_dna_list[i].ntuple(self.virus_dna_list[i]._DNA,N) for i in range(len(self.virus_dna_list))]
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


	def ORFs(self):
		ORF_list=[]
		for item in self.virus_dna_list:
			item.ORF()
			ORF_list.append(item.ORFs)
		return ORF_list

	def ORF_count(self,na_comb):
		count=[]
		N=len(na_comb[0])
		ORFs=self.ORFs()
		if not len(ORFs):
			return []
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
			result.append([[count[i][Tup]] for i in range(len(count))])
		return result
