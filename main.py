# using the dnaAnalyzer class

#to import the class:

from dnaAnalyzerClass import DnaAnalyzer

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

#Print output
dna.reportOutput()

#Write output to given path
dna.writeToFile('EGFPoutput.txt')