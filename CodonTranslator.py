import itertools
import more_itertools
from itertools import zip_longest as izip_longest

class Codons:

    dna_table = {'TTT':'Phe', 'TCT':'Ser', 'TAT':'Tyr', 'TGT':'Cys',
                 'TTC':'Phe', 'TCC':'Ser', 'TAC':'Tyr', 'TGC':'Cys',
                 'TTA':'Leu', 'TCA':'Ser', 'TAA':'Stop', 'TGA':'Stop',
                 'TTG':'Leu', 'TCG':'Ser', 'TAG':'Stop', 'TGG':'Trp',
                 'CTT':'Leu', 'CCT':'Pro', 'CAT':'His', 'CGT':'Arg',
                 'CTC':'Leu', 'CCC':'Pro', 'CAC':'His', 'CGC':'Arg',
                 'CTA':'Leu', 'CCA':'Pro', 'CAA':'Gln', 'CGA':'Arg',
                 'CTG':'Leu', 'CCG':'Pro', 'CAG':'Gln', 'CGG':'Arg',
                 'ATT':'Ile', 'ACT':'Thr', 'AAT':'Asn', 'AGT':'Ser',
                 'ATC':'Ile', 'ACC':'Thr', 'AAC':'Asn', 'AGC':'Ser',
                 'ATA':'Ile', 'ACA':'Thr', 'AAA':'Lys', 'AGA':'Arg',
                 'ATG':'Met', 'ACG':'Thr', 'AAG':'Lys', 'AGG':'Arg',
                 'GTT':'Val', 'GCT':'Ala', 'GAT':'Asp', 'GGT':'Gly',
                 'GTC':'Val', 'GCC':'Ala', 'GAC':'Asp', 'GGC':'Gly',
                 'GTA':'Val', 'GCA':'Ala', 'GAA':'Glu', 'GGA':'Gly',
                 'GTG':'Val', 'GCG':'Ala', 'GAG':'Glu', 'GGG':'Gly'}

    rna_table = {'UUU':'Phe', 'UCU':'Ser', 'UAU':'Tyr', 'UGU':'Cys',
                 'UUC':'Phe', 'UCC':'Ser', 'UAC':'Tyr', 'UGC':'Cys',
                 'UUA':'Leu', 'UCA':'Ser', 'UAA':'Stop', 'UGA':'Stop',
                 'UUG':'Leu', 'UCG':'Ser', 'UAG':'Stop', 'UGG':'Trp',
                 'CUU':'Leu', 'CCU':'Pro', 'CAU':'His', 'CGU':'Arg',
                 'CUC':'Leu', 'CCC':'Pro', 'CAC':'His', 'CGC':'Arg',
                 'CUA':'Leu', 'CCA':'Pro', 'CAA':'Gln', 'CGA':'Arg',
                 'CUG':'Leu', 'CCG':'Pro', 'CAG':'Gln', 'CGG':'Arg',
                 'AUU':'Ile', 'ACU':'Thr', 'AAU':'Asn', 'AGU':'Ser',
                 'AUC':'Ile', 'ACC':'Thr', 'AAC':'Asp', 'AGC':'Ser',
                 'AUA':'Ile', 'ACA':'Thr', 'AAA':'Lys', 'AGA':'Arg',
                 'AUG':'Met', 'ACG':'Thr', 'AAG':'Lys', 'AGG':'Arg',
                 'GUU':'Val', 'GCU':'Ala', 'GAU':'Asp', 'GGU':'Gly',
                 'GUC':'Val', 'GCC':'Ala', 'GAC':'Asp', 'GGC':'Gly',
                 'GUA':'Val', 'GCA':'Ala', 'GAA':'Glu', 'GGA':'Gly',
                 'GUG':'Val', 'GCG':'Ala', 'GAG':'Glu', 'GGG':'Gly'}

    def __init__(self, sequence, seqtype):
        self.seq = sequence.upper()
        self.seqtype = seqtype

    # Groups Codons into tuples
    def grouper(iterable, n, fillvalue=None):
	    args = [iter(iterable)] * n
	    return izip_longest(fillvalue=fillvalue, *args)
     
    # Combines Condon tuples into 3 char strings
    def framer(seq):
	    frames = []
	    for i in seq:
		    x = ''.join(i)
		    frames.append(x)
	    return frames

    #frame_1 = list(grouper(sequence[0:], 3, '-'))
    #frame_2 = list(grouper(sequence[1:], 3, '-'))
    #frame_3 = list(grouper(sequence[2:], 3, '-'))
    #codons_1 = framer(frame_1)
    #codons_2 = framer(frame_2)
    #codons_3 = framer(frame_3)
    codons_1 = framer(list(grouper(sequence[0:].upper(), 3, '-')))
    codons_2 = framer(list(grouper(sequence[1:].upper(), 3, '-')))
    codons_3 = framer(list(grouper(sequence[2:].upper(), 3, '-')))

    def print_frames(self):
        print(self.codons_1, self.codons_2, self.codons_3, sep='\n')

    
    def translate(self, codon_list):
        aminos = []
        for codon in codon_list:
            if self.seqtype is 'RNA':
                aminos.append(self.rna_table.get(codon))
            elif self.seqtype is 'DNA':
                aminos.append(self.dna_table.get(codon))
        return aminos, codon_list


    def print_frame(self, codon_list):
        print(self.translate(codon_list)[0], self.translate(codon_list)[1], sep='\n')


    # Input the reading frame which you wish to print: 1, 2, or 3
    def print_frame2(self, frame):
        if frame == 1:
            print(self.translate(self.codons_1)[0], self.translate(self.codons_1)[1], sep='\n')
        if frame == 2:
            print(self.translate(self.codons_2)[0], self.translate(self.codons_2)[1], sep='\n')
        if frame == 3:
            print(self.translate(self.codons_3)[0], self.translate(self.codons_3)[1], sep='\n')
        
