#!/usr/bin/env python3
import argparse
import sys,os
import subprocess
import shutil
import pandas as pd
from Bio import SeqIO
from Bio.SeqIO import FastaIO

#Rename RSEM and/or fasta 
parser = argparse.ArgumentParser(description = "Script to rename RSEM and fasta contig names to the toxin/nontoxin and count number. Uses the 'Counts' sheet")
parser.add_argument("-r","--rsem",
								type=str,
								help="RSEM excel file, required.")
parser.add_argument("-f","--fasta",
									type=str,
									help="Fasta")
parser.add_argument("-o","--output",
									type=str,
									help="Name of output csv and fasta")
parser.add_argument("-nt","--nontoxins_rename",
								action="store_true",
								default=False,
								help="Set True if nontoxins to be fully renamed, not just Nontoxin-1")
args = parser.parse_args()
rsem=args.rsem
f = args.fasta
nt = args.nontoxins_rename

#which ya doing
if rsem == None:
	print("Missing rsem output")
	quit()
if f == None:
	print("Only renaming rsem")
	
df = pd.read_excel(rsem, sheet_name="Counts")
tc_list = df.toxin_class.unique()

for xin in tc_list:
    n=1
    for index,row in df.iterrows():
        if xin == row["toxin_class"]:
            if row["toxin_class"] == "Nontoxin":
                if row["gene_id"].split('|')[2] != 'sp':
                    df.at[index,"short_id"]=(row["gene_id"].split('|')[2])
                else:
                    df.at[index,"short_id"]=(row["gene_id"].split('|')[4])
            else:
                df.at[index,"short_id"]= xin + "-" + str(n)
                n=n+1

if nt == True:
	for ntx in df.loc[df['toxin_class']== "Nontoxin", 'short_id'].unique():
		n=1
		for index,row in df.iterrows():
			df.at[index,"short_id"] = ntx + "-" + str(n)
			n=n+1

unsorted = list(SeqIO.parse(f, 'fasta'))
sequences = [f for f in sorted(unsorted, key=lambda x : x.name)]
df= df.sort_values(by="gene_id")
df = df.reset_index(drop=True)

for i in range (0,len(sequences)):
    if sequences[i].id != df.at[i,"gene_id"]:
        print("not the same!")
        print(sequences[i].id + "is not" + df.at[i, "gene_id"])
    sequences[i].id=df.at[i,"short_id"]
    sequences[i].description = str("")
fastaname = args.output + ".fasta"
fasta_out=FastaIO.FastaWriter(fastaname,wrap=None)
fasta_out.write_file(sequences)

excelname = args.output + ".xlsx"
writer= pd.ExcelWriter(excelname, engine='xlsxwriter')
df.to_excel(writer, sheet_name="test", index=False)
writer.save()
writer.close()