#dnaAnalyzer.py

class DnaAnalyzer:
	
	def __init__(self):
		
		"""
		
		DnaAnalyzer class will analyze the number of nucleoides, codons, amino acids and GC content of our DNA
		
		"""
		
		self.ps = '%'
		
		self.codonTable = { "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
							"UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
							"UAU":"Y", "UAC":"Y", "UAA":"-", "UAG":"-",
							"UGU":"C", "UGC":"C", "UGA":"-", "UGG":"W",
							"CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
							"CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
							"CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
							"CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
							"AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
							"ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
							"AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
							"AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
							"GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
							"GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
							"GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
							"GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",
							}
		
		#initialize our codon and amino acid dictionaries
		self.codon = {}
		self.aa = {}
		
		self.nuc = { "A":0, "T":0, "G":0, "C":0, "U":0, "N":0 }
		
		#add keys to codon and aa, and init to zero
		for key in self.codonTable:
			self.codon[key] = 0
			self.aa[self.codonTable[key]] = 0
	
		#a list to hold our headers
		self.header = []
	
	def readFasta(self, fp):
		"""
		readFasta reads in the fasta,
		calls parseFasta and then calls analyzeSequence
		to analyze the dna
		"""
		
		for head, seq in self.parseFasta(fp):
			#analyzing the sequence
			self.analyzeSequence(seq)
			#saving the header
			if head == '':
				continue
			else:	
				self.header.append(head)
			
	def parseFasta(self, fp):
		'''
		Reads fasta from a file, takes filepath as argument
		outputs header and sequence
		'''
		
		fh = open(fp, 'r')
		
		header = ''
		sequence = ''
		
		for line in fh:
			
			#if we hit a header
			if line.startswith('>'):
				#return the header, seq
				yield header, sequence
				
				#grab the new header
				header = line.replace('>', '')
				
				#reset our sequence string
				sequence = ''
			
				#don't add the dna this round
				continue
			
			sequence += line.strip()
		
		#return the final header and sequence
		yield header, sequence
	
	def convertDNAtoRNA(self, seq):
		'''
		Converts RNA to DNA, also parses the DNA
		'''
		#make seq uppercase
		seq = seq.upper()
		
		temp = ''
		
		#take out all non dna / rna characters
		for base in seq:
			if base in self.nuc:
				temp += base
		
		seq = temp
		
		#convert dna to rna
		seq = seq.replace('T', 'U')
		
		return(seq)
		
	
	def analyzeSequence(self, seq):
		'''
		analyzeSequence increments the respective dictionaries
		'''
		
		#make the sequence uppercase
		seq = seq.upper()
		
		#increment the bases to our nucleotide dictionary
		for base in seq:
			if base in self.nuc:
				self.nuc[base] += 1
		
		#convert to RNA
		rna = self.convertDNAtoRNA(seq)
		
		#add the codons to the codon dictionary
		for i in range(0, len(rna), 3):
			if rna[i:i+3] in self.codonTable:
				self.codon[rna[i:i+3]] += 1
				self.aa[self.codonTable[rna[i:i+3]]] += 1
				
	def gcContent(self):
		
		gc = ( self.nuc["G"] + self.nuc["C"] ) / ( self.nuc["G"] + self.nuc["C"] + self.nuc["T"] + self.nuc["A"] )
	
		gc = gc * 100
		
		return(gc)
	
	def getMB(self):
		
		totalBase = 0.0
		
		for base in self.nuc:
			totalBase += self.nuc[base]
			
		MB = totalBase / 1000000.00
		
		return(MB)
		
	def totalBases(self):
		
		totalBase = 0.0
		
		for base in self.nuc:
			totalBase += self.nuc[base]
			
		return(totalBase)
	
	def totalAA(self):
		
		totalAA = 0
		
		for aa in self.aa:
			totalAA += self.aa[aa]
			
		return(totalAA)
	
	def reportAA(self):
		
		keys = list(self.aa.keys())
		keys.sort()
		
		print("Amino Acid Composition")
		
		for key in keys:
			print('%s : %3.2f%s %d ' % (key, self.aa[key] / self.totalAA() * 100, self.ps, self.aa[key] ))
		
		print()
	
	def reportCodon(self):
		
		keys = list(self.codon.keys())
		keys.sort()
		
		print("Codon Usage")
		
		for key in keys:
			print("%s : %s %3.2f %s %d " %(key, self.codonTable[key],
			self.codon[key] / self.totalAA() * 100, self.ps,
			self.codon[key] ))
		
		print()
		
	def reportOutput(self):
		
		print()
		print(self.header[0])
		print()
		print("Total Bases = %d " % self.totalBases())
		print()
		print("MB = %.6f " % self.getMB())
		print()
		print("GC content = %.2f%s " % (self.gcContent(), self.ps))
		print()
		
		self.reportAA()
		
		self.reportCodon()
		
	# (Helper) Write the output to a text file
	def writeToFile(self, path):
		f = open(path, 'w')
		
		f.write(self.header[0])
		f.write("\n")
		f.write("Total Bases = %d " % self.totalBases())
		f.write("\n")
		f.write("\n")
		f.write("MB = %.6f " % self.getMB())
		f.write("\n")
		f.write("\n")
		f.write("GC content = %.2f%s " % (self.gcContent(), self.ps))
		
		#Write AA to file
		keys = list(self.aa.keys())
		keys.sort()
		
		f.write("\n\nAmino Acid Composition\n")
		
		for key in keys:
			f.write('%s : %3.2f%s %d \n' % (key, self.aa[key] / self.totalAA() * 100, self.ps, self.aa[key] ))
		
		print("Done writing AA to file!")
		
		#Write Codon to file
		keys = list(self.codon.keys())
		keys.sort()
		
		f.write("\nCodon Usage\n")
		
		for key in keys:
			f.write("%s : %s %3.2f %s %d \n" %(key, self.codonTable[key],
			self.codon[key] / self.totalAA() * 100, self.ps,
			self.codon[key] ))
		
		f.write("\n")
		
		print("Done writing Codon to file!")

		print("Done writing to file!")
		f.close

		

	
#Create a DnaAnalyzer Object
dna = DnaAnalyzer()

#Read the egfp.fasta file
'''
Enhanced green fluorescent protein
is a protein composed of 238 amino acid residues (26.9kDa)
that exhibits bright green fluorescence when exposed to
light in the blue to ultraviolet range.
'''
dna.readFasta('egfp.fasta')
dna.reportOutput()














'''	
# (Helper) Write the output to a text file
def writeToFile(aDict, path):
	f = open(path, 'w')
	for key in sorted(aDict):
		f.write("%s: %s" % (key, aDict[key]))
		f.write("\n")
	print("Done writing to file!")
	f.close	
'''

'''
#Write each output to their corresponding file
writeToFile(dna.codon, "codon.txt")
writeToFile(dna.aa, "aa.txt")
writeToFile(dna.nuc, "nuc.txt")
'''

