from collections import Counter

class virus:
	def __init__(self,virus_dna):
		if(type(virus_dna)==dict):
			self._DNA=virus_dna["RNA"]
			self._NAME=virus_dna["NAME"]
		
	def ntuple(self,N):
		Split_DNA=[self._DNA[i:i+N] for i in range(0,len(self._DNA),N)]
		Tuplecomb = dict(Counter(Split_DNA))
		return Tuplecomb

class virus_species:
	def __init__(self,virus_dna):
		if(type(virus_dna)==list):
			self.virus_dna_list=[]
			for element in virus_dna:
				self.virus_dna_list.append(virus(element))
	
	def ntuple(self,N):
		tuple_list = [self.virus_dna_list[i].ntuple(N) for i in range(len(self.virus_dna_list))]
		return tuple_list
	
	def na_count(self,na_comb,na_comb1):
		tuplelist=self.ntuple(len(na_comb))
		count_0 = [tuplelist[i][na_comb] for i in range(len(tuplelist))]
		count_1 = [tuplelist[i][na_comb1] for i in range(len(tuplelist))]	
		return count_0,count_1
