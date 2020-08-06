#!/usr/bin/python
#coding: utf-8 

#######################################################################
#translates a DNA sequence, in FASTA format, into six protein sequences
#based on the six open reading frames
#######################################################################

import re

def open_file(fasta_to_open):
    """
    opens fasta file and creates dictionary of chromosomal sequences
    ensures all values are lower for conversion function
    creates a dictionary with the chromosome number as the key and its sequence as the value
    """
    fasta = open(fasta_to_open, "r")
    fasta = fasta.read().split("\n\n")
    fasta = [fasta_line.split("\n",1) for fasta_line in fasta]
    for chr_line in fasta:
        chr_seq = chr_line[1]
        chr_seq = chr_seq.replace("\n","")
        chr_seq = chr_seq.upper()
    return chr_seq

def rev_comp_func(seq):
    """ generates the reverse complement of a DNA sequence"""
    nt = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nt[n] for n in reversed(seq))
    
def codon_func(seq):
    """ splits each mRNA sequence by 3 letters and removes any string 
    not divisible by 3 to ensure each codon is 3 letters long"""
    codons_string = " ".join(seq[i:i+3] for i in range(0, len(seq), 3))
    codons_list = codons_string.split(" ")
    codons = []
    for a in codons_list:
            if len(a) % 3 == 0:
                codons.append(a)
    return codons

def protein_trans_func(codons):
    """ takes a list of codons, translates each codon to an amino acid/
    stop codon, merges the product back together, and only keeps translated 
    products that begin with Methionine and end with a stop codon"""
    prot_list = []
    alanine = [re.sub("(GC)[AGCU]", "A", i) for i in codons]
    arginine = [re.sub("(CG)[AGCU]","R", i) for i in alanine]
    arginine = [re.sub("(AG)[AG]", "R", i) for i in arginine]
    asparagine = [re.sub("(AA)[UC]", "N", i) for i in arginine]
    aspartic_acid = [re.sub("(GA)[UC]", "D", i) for i in asparagine]
    cysteine = [re.sub("(UG)[UC]", "C", i) for i in aspartic_acid]
    glutamine = [re.sub("(CA)[AG]", "Q", i) for i in cysteine]
    glutamic_acid = [re.sub("(GA)[AG]", "E", i) for i in glutamine]
    glycine = [re.sub("(GG)[AGCU]", "G", i) for i in glutamic_acid]
    histidine = [re.sub("(CA)[UC]", "H", i) for i in glycine]
    isoleucine = [re.sub("(AU)[UCA]", "I", i) for i in histidine]
    leucine = [re.sub("(CU)[AGCU]", "L", i) for i in isoleucine]
    leucine = [re.sub("(UU)[AG]", "L", i) for i in leucine]
    lysine = [re.sub("(AA)[AG]", "K", i) for i in leucine]
    methionine = [re.sub("(AUG)", "M", i) for i in lysine]
    phenylalanine = [re.sub("(UU)[UC]", "F", i) for i in methionine]
    proline = [re.sub("(CC)[AGCU]", "P", i) for i in phenylalanine]
    serine = [re.sub("(UC)[AGCU]", "S", i) for i in proline]
    serine = [re.sub("(AG)[UC]", "S", i) for i in serine]
    threonine = [re.sub("(AC)[AGCU]", "T", i) for i in serine]
    tryptophan = [re.sub("(UGG)", "W", i) for i in threonine]
    tyrosine = [re.sub("(UA)[UC]", "Y", i) for i in tryptophan]
    valine = [re.sub("(GU)[AGCU]", "V", i) for i in tyrosine]
    stop_codon = [re.sub("(UA)[GA]", "-", i) for i in valine]
    stop_codon = [re.sub("(UGA)", "-", i) for i in stop_codon]
    clean_seq = "".join(stop_codon)
    split_by_end = clean_seq.replace("-","- ").split(" ")
    split_by_start = [aa.replace("M"," M",1).split(" ") for aa in split_by_end]
    proteins = [seq for sublist in split_by_start for seq in sublist]
    [prot_list.append(i) for i in proteins if i.startswith("M") and i.endswith("-")]
    return prot_list

# generate the first three frames for each DNA sequence:
dna_seqs1 = open_file(path)
dna_seqs2 = dna_seqs1[1:]
dna_seqs3 = dna_seqs2[1:]

# generate the reverse complements for each of those frames:
dna_rev_comps1 = rev_comp_func(dna_seqs1)
dna_rev_comps2 = dna_rev_comps1[1:]
dna_rev_comps3 = dna_rev_comps2[1:]

# let's transcribe these sequences to RNA:
rna_seqs1 = re.sub("T","U",dna_seqs1)
rna_seqs2 = re.sub("T","U",dna_seqs2)
rna_seqs3 = re.sub("T","U",dna_seqs3)
rna_rev_comps1 = re.sub("T","U",dna_rev_comps1)
rna_rev_comps2 = re.sub("T","U",dna_rev_comps2)
rna_rev_comps3 = re.sub("T","U",dna_rev_comps3)

# now that we have our mRNA we can separate each sequence into codons:
codon_seqs1 = codon_func(rna_seqs1)
codon_seqs2 = codon_func(rna_seqs2)
codon_seqs3 = codon_func(rna_seqs3)
codon_rev_comps1 = codon_func(rna_rev_comps1)
codon_rev_comps2 = codon_func(rna_rev_comps2)
codon_rev_comps3 = codon_func(rna_rev_comps3)

# last step! now each codon needs to be translated and the translated product
# needs to be trimmed so that it starts with a start codon (Methionine) and 
# ends with a stop codon not to be translated (UAG, UAA, and UGA) [2]:
protein_1 = protein_trans_func(codon_seqs1)
protein_2 = protein_trans_func(codon_seqs2)
protein_3 = protein_trans_func(codon_seqs3)
protein_4 = protein_trans_func(codon_rev_comps1)
protein_5 = protein_trans_func(codon_rev_comps2)
protein_6 = protein_trans_func(codon_rev_comps3)

# let's print out the starting sequence and the resulting amino acid sequence:
print(dna_seqs1)
print(protein_1)
print(dna_seqs2)
print(protein_2)
print(dna_seqs3)
print(protein_3)
print(dna_rev_comps1)
print(protein_4)
print(dna_rev_comps2)
print(protein_5)
print(dna_rev_comps3)
print(protein_6)

# References:

# [1] Clancy, S. (2008) RNA splicing: introns, exons and spliceosome.
#     Nature Education 1(1):31. https://www.nature.com/scitable/topicpage/rna-splicing-introns-exons-and-spliceosome-12375/
# [2] Genetic Code. Retrieved on 04 August 2020 from https://en.wikipedia.org/wiki/Genetic_code
