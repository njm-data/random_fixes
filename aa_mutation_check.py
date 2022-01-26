#!/usr/bin/env python3
import argparse
import sys,os
import subprocess
import shutil
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqIO import FastaIO

#Tells you what amino acid family mutations fall into
parser = argparse.ArgumentParser(description = "Tells you what amino acid family mutations fall into")
parser.add_argument("-f","--fasta",
								type=str,
								help="Translated and aligned fasta")
parser.add_argument("-s","--single",
								action="store_true",
								default=False,
								help="Will compare all sequences to the first sequence in the alignment as the reference")
parser.add_argument("-p","--pairwise",
								action="store_true",
								default=False,
								help="Will compare all sequences to each other in the alignment")

args = parser.parse_args()
f = args.fasta
single = args.single
pairwise = args.pairwise

if f == None:
	print("Missing fasta")
	quit()
if single and pairwise == None:
	print("Must designate one-way or pairwise comparisons")
	quit()


aa_list = [{'positve_aa': ['R', 'H', 'K']}, #Arginie, Histidine, Lysine
           {'negative_aa': ['D', 'E']}, # Aspartic Acid, Glutamic Acid
           {'polar_un': ['S', 'T', 'N', 'Q']}, #Serine, Threonine, Asparagine, Glutamine
           {'hydrophobic' : ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']}, #Alanine, Valine, Isoleucine, Leucine, Methionine, 
                                                        #Phenylalanine, Tyrosine, Trytophan 
           {'special_aa' : ['C', 'G', 'P', 'U']}, #Cysteine, Glycine, Proline, Selenocysteine
           {'indel' : '-'}]
           
sequences = list(SeqIO.parse(f,'fasta'))

if single:
    reference = sequences[0].seq
    reflen= len(reference) -1
    for x in range (1,len(sequences)): #go through fasta
        n = 0
        sametally = 0
        ref_res = []
        new_res = []
        comparison = ("Comparing " + sequences[0].id + " and " + sequences[x].id)
        print(comparison)
        for y in sequences[x].seq: #go through individual residues
            if n > reflen: #if reference is shorter than other sequences, stops here.
                print("reference is too short, something is misaligned")
                break
            
            if y != reference[n]: #if there's a mutation
                for dict_item in aa_list: #go through amino acid list and catalog the mutation types
                    for key in dict_item:
                        for residue in dict_item[key]:
                            if residue == reference[n]:
                                a = key
                                ref_res.append(key)
                            if residue == y:
                                b = key
                                new_res.append(key)
                change = str("Went from " + reference[n] + " to a " + y + " at position " + str(n) + "\n" +
                            a + " to " + b)
                print(change)
            else:
                sametally = sametally + 1
            n = n + 1
        if ref_res:
            for dict_item in aa_list:
                for key in dict_item:
                    print(sequences[0].id + " (reference) had " + str(ref_res.count(key)) + " " + key + " residues")
                    print(sequences[x].id + " had "+ str(new_res.count(key)) + " " + key + " residues")
            print(str(sametally) + " residues were the same.\n")    
        else:
            print("Identical\n")


if pairwise:
	for r in range(0,len(sequences)):
		reference = sequences[r].seq
		reflen= len(reference) -1
		for x in reversed(range(0,len(sequences))): #go through fasta
			if r == x:
				break
			n = 0
			sametally = 0
			ref_res = []
			new_res = []
			comparison = ("Comparing " + sequences[r].id + " and " + sequences[x].id)
			print(comparison)
			for y in sequences[x].seq: #go through individual residues
				if n > reflen: #if reference is shorter than other sequences, stops here.
					print("reference is too short, something is misaligned")
					break
			
				if y != reference[n]: #if there's a mutation
					for dict_item in aa_list: #go through amino acid list and catalog the mutation types
						for key in dict_item:
							for residue in dict_item[key]:
								if residue == reference[n]:
									a = key
									ref_res.append(key)
								if residue == y:
									b = key
									new_res.append(key)
					change = str("Went from " + reference[n] + " to a " + y + " at position " + str(n+1) + "\n" +
								a + " to " + b)
					print(change)
				else:
					sametally = sametally + 1
				n = n + 1
			if ref_res:
				for dict_item in aa_list:
					for key in dict_item:
						print(sequences[r].id + " (reference) had " + str(ref_res.count(key)) + " " + key + " residues")
						print(sequences[x].id + " had "+ str(new_res.count(key)) + " " + key + " residues")
				print(str(sametally) + " residues were the same.\n")    
			else:
				print("Identical\n")