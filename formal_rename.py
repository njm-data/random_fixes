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
								help="RSEM excel file")
parser.add_argument("-f","--fasta",
									type=str,
									help="Fasta")
parser.add_argument("-o","--output",
									type=str,
									help="Name of output csv and fasta")
parser.add_argument("-nt","--nontoxins_rename",
								action="store_true",
								default=False,
								help="Have if nontoxins to be fully renamed, not just Nontoxin-1")
args = parser.parse_args()
rsem=args.rsem
f = args.fasta
nt = args.nontoxins_rename

#what's missing
if rsem == None:
	print("Missing rsem output")
	quit()
if f == None:
	print("missing fasta")
	
df = pd.read_excel(rsem, sheet_name="Counts")
sequences = list(SeqIO.parse(f, 'fasta'))

if len(df) != len(sequences):
	print("RSEM and FASTA are different lengths \n RSEM is " + str(len(df)) + "\n FASTA is " + str(len(sequences)))
	quit()

df["average"]=df.mean(axis=1)
df = df.sort_values(by=["average"], ascending=False).reset_index(drop=True)
tc_list = df.toxin_class.unique()

for xin in tc_list:
    n=1
    for index,row in df.iterrows():
        if xin == row["toxin_class"]:
            	df.at[index,"short_id"]= xin + "-" + str(n)
            	n=n+1

if nt == True:
	for index,row in df.iterrows():
		if row["toxin_class"] == "Nontoxin":
			if row["gene_id"].split('|')[2] != "sp":
				df.at[index,"short_id"]=(row["gene_id"].split('|')[2])
			else:
				df.at[index,"short_id"]=(row["gene_id"].split('|')[4])
	for ntx in df.loc[df['toxin_class']== "Nontoxin", 'short_id'].unique():
		n=1
		for index,row in df.iterrows():
			if ntx == row["short_id"]:
				df.at[index,"short_id"] = ntx + "-" + str(n)
				n=n+1

sequences = list(SeqIO.parse(f, 'fasta'))
sortedList = [f for f in sorted(sequences, key=lambda x : x.name)]
df1 = df.sort_values(by=["gene_id"]).reset_index(drop=True)
for i in range (0,len(sortedList)):
    sortedList[i].id=df1.at[i,"short_id"]
    sortedList[i].description = str("")
fastaname = args.output + ".fasta"
fasta_out=FastaIO.FastaWriter(fastaname,wrap=None)
fasta_out.write_file(sortedList)

excelname = args.output + ".xlsx"
writer= pd.ExcelWriter(excelname, engine='xlsxwriter')
df.to_excel(writer, sheet_name="test", index=False)
writer.save()
writer.close()
