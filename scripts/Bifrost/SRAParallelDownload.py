from Bio import Entrez
import csv
import pandas as pd
import numpy as np
import warnings
from multiprocessing import Pool
import concurrent.futures
import multiprocessing as mp
import os
import shutil
import subprocess

warnings.simplefilter(action='ignore', category=FutureWarning)
Entrez.email = "rugarem@liverpool.ac.uk"

#load the original csv file to be downsampled
df = pd.read_csv('/data/sra/hydra2/vus/toolkits/sra/download_runs.csv', index_col=False)
np.random.seed(100)

# Create a new dataframe to store the downsampled runs
selected_df = pd.DataFrame(columns=df.columns)
# Filter out runs larger than 10GB
filtered_df = df[df['Size (MB)'] <= 10240]


# Loop through each unique tissue and randomly select run IDs until the total size is up to 20GB
for tissue in df['Tissue'].unique():
    print (tissue)
    tissue_df = filtered_df[filtered_df['Tissue'] == tissue]
    tissue_df = tissue_df.sample(frac=1, random_state=100)
    total_size = 0
    for _, row in tissue_df.iterrows():
        if total_size + row['Size (MB)'] <= 20480:
            #print(row)
            selected_df = selected_df.append(row)
            total_size += row['Size (MB)']
        else:
            break
selected_df.to_csv("downsampled_runs.csv", index=False)

#make a dict of the SRR IDs and the corresponding tissue code
dictionary = {}
with open("downsampled_runs.csv", newline='') as csvfile:
    reader = csv.reader(csvfile)
    #skip the header
    header = next(reader)
    dictionary = {row[0]: row[4] for row in reader} 

#print(dictionary)	
key_list = list(dictionary.keys())

new_names = []

#in parallel batches of n, download the fastq files and rename them ... 
def process_srr(srr):
    subprocess.call (['/data/sra/hydra2/vus/toolkits/bifrost/sratoolkit.3.0.5-ubuntu64/bin/prefetch', srr])
    subprocess.call (['/data/sra/hydra2/vus/toolkits/bifrost/sratoolkit.3.0.5-ubuntu64/bin/fasterq-dump','-v', "--concatenate-reads", '--threads', '100', srr])
    
    #remove sra directory
    tmp_path = "/home/rmaruzani/storage/hydra2/vus/toolkits/bifrost/parallel"
    for d in os.listdir(tmp_path):
        if d.endswith(srr):
            print("removing prefetch directory " + d, flush=True)
            shutil.rmtree(os.path.join(tmp_path, d))
            
    current_name = srr + ".fastq"
    new_name = dictionary[srr]+"_"+current_name
    os.rename(current_name, new_name)
    print ("renaming fastq " + current_name + " to " + new_name + "\n", flush=True)
    return new_name

#run parallel processes
if __name__ == '__main__':
    srrs = list(dictionary.keys())
    pool = Pool(processes=50)  # create a pool of 10 worker processes
    new_names = pool.map(process_srr, srrs)  #process 10 SRRs at a time in parallel
    pool.close()
    pool.join() 

#write the new names to a file as input for Bifrost    
with open('sra.txt', 'w') as f:
    for srr in new_names:
        f.write(srr+'\n')

#run Bifrost build on the downloaded fastq files
with open('bifrost_output.log', 'a') as log_file:
    subprocess.call(["Bifrost", "build", "-v", "-t", "150", "-k", "31", "-c", "-i", "-d", "-s", "sra.txt", "-r", "/data/sra/hydra2/vus/toolkits/bifrost/hg38.fa", "-o", "metasra"], stdout=log_file, stderr=log_file)

