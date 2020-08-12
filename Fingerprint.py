#!/usr/bin/python3
from virus import virus, virus_species, readfile
import pickle as pk
import sys

PROT={'TTT':'P','TTC':'P','TTA':'L','TTG':'L','CTT':'L','CTC':'L',
'CTA':'L','CTG':'L','ATT':'Is','ATC':'Is','ATA':'Is','ATG':'M',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','TCT':'S','TCC':'S',
'TCA':'S','TCG':'S','CCT':'Pr','CCC':'Pr','CCA':'Pr','CCG':'Pr',
'ACT':'Th','ACC':'Th','ACA':'Th','ACG':'Th','GCT':'Al',
'GCC':'Al','GCA':'Al','GCG':'Al','TAT':'Ty','TAC':'Ty',
'TAA':'STOP','TAG':'STOP','CAT':'H','CAC':'H','CAA':'Gl',
'CAG':'Gl','AAT':'As','AAC':'As','AAA':'Ly','AAG':'Ly',
'GAT':'Asp','GAC':'Asp','GAA':'Glu','GAG':'Glu','TGT':'C',
'TGC':'C','TGA':'STOP','TGG':'Try','CGT':'Arg','CGC':'Arg',
'CGA':'Arg','CGG':'Arg','AGT':'S','AGC':'S','AGA':'Arg',
'AGG':'Arg','GGT':'Gly','GGC':'Gly','GGA':'Gly','GGG':'Gly'}

PROT_code={'P':'TTY','L':['TTR','CTN'],'Is':'ATH','M':'ATG','V':'GTN','S':['TCN','AGY'],'Pr':'CCN','Th':'ACN','Al':'GCN','Ty':'TAY','STOP':['TAR','TGA'],'H':'CTY','Gl':'CAR','As':'AAY','Ly':'AAR','Asp':'GAN','Glu':'GAR','C':'TGY','Try':'TGG','Arg':['CGN','AGR'],'Gly':'GGN'}

resolution = {'A':'A','C':'C','G':'G','T':'T','R':['A','G'],'Y':['C','T'],'W':['A','T'],'S':['G','C'],'M':['A','C'],'K':['G','T'],'H':['A','C','T'],'B':['G','C','T'],'V':['G','A','C'],'D':['G','A','T'],'N':['A','C','G','T']}

class MissingFingerprint(Exception):
     def __str__(self):
        return "Fingerprint consists of only wildcards or no fingerprint was found"

class ShortVirus(Exception):
     def __str__(self):
        return "Virus DNA sequence is to short"


def ORF_read(filename):
	""" Read file and return ORF
		
	A function to read in files with fasta format and return ORF sequences
	in an array.

	args:
		filename (str) -- Name of the fasta file
	
	returns:
		Virus_ORFs (list) -- A list with all virus ORFs from the DNA sequences contained in the file
	"""
	Virus_DNA=readfile(filename)
	Virus_DNA=virus_species(Virus_DNA)
	Virus_ORFs=Virus_DNA.ORFs()
	return Virus_ORFs

def ORF_scan(filename):
	""" Search ORF for Fingerprint

	A Function that searches for Nucleobases that are situated at the same position for all sequences.
	The translation of the sequence is used aswell to determine if some aminoacids occur at the same positions
	to include the corresponding nucleobase triplet as well. It stores the Fingerprint in a pickle file.

	args:
		filename (str) -- The name of the fasta file containing the virus DNA

	returns:
	"""
	ORFs=ORF_read(filename)
	lengths=[]
	for item in ORFs:
		for item1 in item:
			lengths.append(len(item1['ORF']))
	
	minimum = min(lengths)
		

	Fingerprint=[]
	counter=0
	for i in range(0,minimum,3):
		if(ORFs[0]==[]):
			continue
		Fingerprint.append(ORFs[0][0]['ORF'][i])
		Fingerprint.append(ORFs[0][0]['ORF'][i+1])
		Fingerprint.append(ORFs[0][0]['ORF'][i+2])
		for j in range(0,len(ORFs),len(ORFs)//150):
			if(ORFs[j]==[]):
				continue
			if(ORFs[j][0]['ORF'][i]!='A' and ORFs[j][0]['ORF'][i]!='C' and ORFs[j][0]['ORF'][i]!='T' and ORFs[j][0]['ORF'][i]!='G'):
				continue
			if(ORFs[j][0]['ORF'][i+1]!='A' and ORFs[j][0]['ORF'][i+1]!='C' and ORFs[j][0]['ORF'][i+1]!='T' and ORFs[j][0]['ORF'][i+1]!='G'):
				continue
			if(ORFs[j][0]['ORF'][i+2]!='A' and ORFs[j][0]['ORF'][i+2]!='C' and ORFs[j][0]['ORF'][i+2]!='T' and ORFs[j][0]['ORF'][i+2]!='G'):
				continue
			if(PROT[ORFs[0][0]['ORF'][i:i+3]]!=PROT[ORFs[j][0]['ORF'][i:i+3]]):
				Fingerprint[i]='N'
				Fingerprint[i+1]='N'
				Fingerprint[i+2]='N'
				counter+=3
				break
	
			if(Fingerprint[i]!=ORFs[j][0]['ORF'][i]):
				if(type(PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]])==list):
					Fingerprint[i]=[PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][0][0],PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][1][0]]
				else:
					Fingerprint[i]=PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][0]
			if(Fingerprint[i+1]!=ORFs[j][0]['ORF'][i+1]):
				if(type(PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]])==list):
					if(Fingerprint[i]==[PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][0][0],PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][1][0]]):
						Fingerprint[i+1]=[PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][0][1],PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][1][1]]
					if(Fingerprint[i]==PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][0][0]):
						Fingerprint[i+1]=PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][0][1]
					if(Fingerprint[i]==PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][1][0]):
						Fingerprint[i+1]=PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][1][1]
				else:
					Fingerprint[i+1]=PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][1]
			if(Fingerprint[i+2]!=ORFs[j][0]['ORF'][i+2]):
				if(type(PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]])==list):
					if(Fingerprint[i]==[PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][0][0],PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][1][0]]):
						Fingerprint[i+2]=[PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][0][2],PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][1][2]]
					if(Fingerprint[i]==PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][0][0]):
						Fingerprint[i+2]=PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][0][2]
					if(Fingerprint[i]==PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][1][0]):
						Fingerprint[i+2]=PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][1][2]
				else:
					Fingerprint[i+2]=PROT_code[PROT[ORFs[0][0]['ORF'][i:i+3]]][2]
	if(not len(Fingerprint)):
		raise MissingFingerprint
	if(counter==len(Fingerprint)-3):
		raise MissingFingerprint
	pk.dump(Fingerprint,open("Fingerprint",'wb'))




def test_identity(filename):
	"""
	Tests if sequences contain Fingerprint

	A function that tests if sequences from a fasta file contain the fingerprint stored
	in the Fingerprint pickle file. It will test all positions in the fingerprint
	for identity. Only if all positions are in accordance with the fingerprint the
	sequence is categorized into the category of virus the fingerprint stems from.
	The function will print the percentage of sequences containing the fingerprint
	as a result.

	args:
		filename (str) -- Name of the fasta file containing the virus dna sequences to test
	
	returns:
		
	"""
	ORFs=ORF_read(filename)
	Fingerprint=pk.load(open("Fingerprint",'rb'))
	lengths=[]
	for item in ORFs:
		for item1 in item:
			lengths.append(len(item1['ORF']))
	lengths.append(len(Fingerprint))
	minimum = min(lengths)
	maximum = max(lengths)
	if(maximum<len(Fingerprint)):
		raise ShortVirus
	
	fits=[]
	for i in range(len(ORFs)):
		fits.append(1)
		if(ORFs[i] == []):
			fits[i]=0
			continue
		for j in range(minimum):
			if(ORFs[i][0]['ORF'][j] in ['R','Y','W','S','M','K','H','B','V','D','N']):
					#TODO Skip mechanic
					continue
			if(type(Fingerprint[j]) is list):
				if(ORFs[i][0]['ORF'][j] not in resolution[Fingerprint[j][0]] and ORFs[i][0]['ORF'][j] not in resolution[Fingerprint[j][1]]):
					fits[i]=0
					continue
			else:	
				if(ORFs[i][0]['ORF'][j] not in resolution[Fingerprint[j]]):
					fits[i]=0
					continue
	score=0
	for item in fits:
		score+=item
	print(score/len(fits))

if __name__ == '__main__':
	if(len(sys.argv)==1):
		raise ValueError("Please supply a fasta file containing the virus DNA for testing as command line argument")
	test_identity(sys.argv[1])

		
		

