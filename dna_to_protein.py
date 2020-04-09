#!/usr/bin/python
#coding: utf-8 

#######################################################################
#translates a DNA sequence, in FASTA format, into six protein sequences
#based on the six open reading frames
#######################################################################
import re

chr_dict = {}
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
        marker = chr_line[0]
        chr_seq = chr_line[1]
        chr_seq = chr_seq.replace("\n","")
        chr_seq = chr_seq.lower()
        chr_dict[marker] = chr_seq
    return chr_dict


comp_dict = {}
def make_complement():
    """
    makes the complement DNA strand
    converts a to t and g to c
    """
    for key, value in chr_dict.items():
        a_to_t = value.replace("a","T")
        t_to_a = a_to_t.replace("t","A")
        g_to_c = t_to_a.replace("g","C")
        c_to_g = g_to_c.replace("c","G")
        comp_val = c_to_g.lower()
        comp_dict[key] = comp_val
        

def reading_frames(original_dict):
    """
    creates three different reading frames for chromosomal sequences
    creates three different dictionaries for those reading frames
    """
    dict1 = {}
    dict2 = {}
    dict3 = {}
    for key, value in original_dict.items():
        frame_1 = value
        frame_2 = value[1:]
        frame_3 = value[2:]
        dict1[key] = frame_1
        dict2[key] = frame_2
        dict3[key] =frame_3
    return dict1,dict2,dict3

def codon_gen(frame_dict):
    """
    generates a list of codons from a particular reading frame     
    inserts space after each codon and splits codons into a list of codons
    only appends codons to list of codons
    turns list of codons into codon string
    creates a codon dictionary
    """
    codon_dict = {}
    for key, value in frame_dict.items():
        val_list = []
        val = " ".join(value[i:i+3] for i in range(0, len(value), 3))
        val_list = val.split(" ")
        codon_list = []
        for a in val_list:
            if len(a) % 3 == 0:
                codon_list.append(a)
        codon_dict[key] = codon_list
    return codon_dict


def protein_gen(codon_dict):
    """
    converts codons to proteins
    only takes proteins that start with M and stop with a stop codon
    protein conversions based on genetic table from wiki
    """
    protein_dict = {}
    for key, value in codon_dict.items():
        prot_list = []
        alanine = [re.sub("(gc)[agct]", "A", i) for i in value]
        arginine = [re.sub("(cg)[agct]","R", i) for i in alanine]
        arginine = [re.sub("(ag)[ag]", "R", i) for i in arginine]
        asparagine = [re.sub("(aa)[tc]", "N", i) for i in arginine]
        aspartic_acid = [re.sub("(ga)[tc]", "D", i) for i in asparagine]
        cysteine = [re.sub("(tg)[tc]", "C", i) for i in aspartic_acid]
        glutamine = [re.sub("(ca)[ag]", "Q", i) for i in cysteine]
        glutamic_acid = [re.sub("(ga)[ag]", "E", i) for i in glutamine]
        glycine = [re.sub("(gg)[agct]", "G", i) for i in glutamic_acid]
        histidine = [re.sub("(ca)[tc]", "H", i) for i in glycine]
        isoleucine = [re.sub("(at)[tca]", "I", i) for i in histidine]
        leucine = [re.sub("(ct)[actg]", "L", i) for i in isoleucine]
        leucine = [re.sub("(tt)[ag]", "L", i) for i in leucine]
        lysine = [re.sub("(aa)[ag]", "K", i) for i in leucine]
        methionine = [re.sub("(atg)", "M", i) for i in lysine]
        phenylalanine = [re.sub("(tt)[tc]", "F", i) for i in methionine]
        proline = [re.sub("(cc)[tcag]", "P", i) for i in phenylalanine]
        serine = [re.sub("(tc)[tcag]", "S", i) for i in proline]
        serine = [re.sub("(ag)[tc]", "S", i) for i in serine]
        threonine = [re.sub("(ac)[tcag]", "T", i) for i in serine]
        tryptophan = [re.sub("(tgg)", "W", i) for i in threonine]
        tyrosine = [re.sub("(ta)[tc]", "Y", i) for i in tryptophan]
        valine = [re.sub("(gt)[tcag]", "V", i) for i in tyrosine]
        stop_codon = [re.sub("(ta)[ga]", "*", i) for i in valine]
        stop_codon = [re.sub("(tga)", "*", i) for i in stop_codon]
        clean_seq = "".join(stop_codon)
        start = clean_seq.replace("M"," M").replace("*","* ")
        split = start.split(" ")
        [prot_list.append(i) for i in split if i.startswith("M") and i.endswith("*")]
        protein_dict[key] = prot_list
    return protein_dict

#opens the fasta file
open_file(r"./*.fasta")

#creates three reading frames for the first DNA strand 
norm_frame_one = reading_frames(chr_dict)[0]
norm_frame_two = reading_frames(chr_dict)[1]
norm_frame_three = reading_frames(chr_dict)[2]

#makes the complement DNA strand
make_complement()

#makes the three reading frames for the complement DNA strand
comp_frame_one = reading_frames(chr_dict)[0]
comp_frame_two = reading_frames(chr_dict)[1]
comp_frame_three = reading_frames(chr_dict)[2]

#makes the codon dictionary for the first DNA strand's reading frames
#makes the codon dictionary for the complement DNA strand's reading frames
norm_codon_one = codon_gen(norm_frame_one)
norm_codon_two = codon_gen(norm_frame_two)
norm_codon_three = codon_gen(norm_frame_three)
comp_codon_one = codon_gen(norm_frame_one)
comp_codon_two = codon_gen(norm_frame_two)
comp_codon_three = codon_gen(norm_frame_three)

#makes the orf dictionary for the first DNA strand's reading frames
#makes the orf dictionary for the complement DNA strand's reading frames
norm_orf_one = protein_gen(comp_codon_one)
norm_orf_two = protein_gen(comp_codon_two)
norm_orf_three = protein_gen(comp_codon_three)
comp_orf_one = protein_gen(comp_codon_one)
comp_orf_two = protein_gen(comp_codon_two)
comp_orf_three = protein_gen(comp_codon_three)

print(norm_orf_one)
print(norm_orf_two)
print(norm_orf_three)
print(comp_orf_one)
print(comp_orf_two)
print(comp_orf_three)

    
