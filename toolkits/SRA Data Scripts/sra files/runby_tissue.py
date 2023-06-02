from Bio import Entrez
import subprocess
import os
import csv
import pandas as pd

Entrez.email = "rugarem@liverpool.ac.uk"



#make a dict of the SRR IDs and the corresponding tissue code
dictionary = {}
with open("batch_tissue_runs.csv", newline='') as csvfile:
    reader = csv.reader(csvfile)
    header = next(reader)
    dictionary = {row[0]: row[3] for row in reader} 
	
#Rename the first SRA file, and begin the Bifrost build command, start downloading the next SRA file
for srr in dictionary.keys():
    subprocess.Popen(["./sratoolkit.3.0.2-ubuntu64/bin/fasterq-dump", "-v", "--concatenate-reads", srr], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with open('sra.txt', 'w') as f:
        f.write(srr+".fastq"+'\n')

#subprocess.call(["Bifrost", "build", "-v", "-t", "150", "-k", "31", "-c", "-i", "-d", "-s", "sra.txt", "-r", "/home/rmaruzani/storage/vus/toolkits/bifrost/hg38.fa", "-o", first], stdout=log_file, stderr=log_file)
    