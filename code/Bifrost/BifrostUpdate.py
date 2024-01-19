from Bio import Entrez
import subprocess
import os
import csv
Entrez.email = "rugarem@liverpool.ac.uk"

#make a dict of the SRR IDs and the corresponding tissue code
dictionary = {}
with open('/home/rmaruzani/storage/vus/toolkits/bifrost/batch_runs.csv', newline='') as csvfile:
    reader = csv.reader(csvfile)
    header = next(reader)
    dictionary = {row[0]: row[4] for row in reader} 
	
os.chdir('/home/rmaruzani/storage/vus/toolkits/bifrost')
cwd = os.listdir("/home/rmaruzani/storage/vus/toolkits/bifrost")
first = next(iter(dictionary))
key_list = list(dictionary.keys())

# Download the first SRA file using, rename and build bifrost index
with open('bifrost_output.log', 'a') as log_file:
    print("downloading first fastq  " + first , flush=True)
    subprocess.call(["./sratoolkit.3.0.2-ubuntu64/bin/fasterq-dump", "-v", "--concatenate-reads", first], stdout=log_file, stderr=log_file)
    
    fq_name = first + ".fastq"
    n_name = dictionary[first]+"_"+fq_name
    print ("renaming fastq " + fq_name + " to " + n_name + "\n", flush=True)
    os.rename(fq_name, n_name)
    
    print("building bifrost index with the fastq " + n_name , flush=True)
    subprocess.call(["Bifrost", "build", "-v", "-t", "60", "-k", "31", "-c", "-i", "-d", "-s", n_name, "-r", "/home/rmaruzani/storage/vus/toolkits/bifrost/hg38.fa", "-o", first], stdout=log_file, stderr=log_file)
    print("removing fastq from directory ...", flush=True)
    subprocess.call(["rm", n_name])
    print("removing first accession from dictionary ...", flush=True)
    del dictionary[first]

    
#Rename the first SRA file, and begin the Bifrost build command, start downloading the next SRA file
for srr in dictionary.keys():
    try:
        next_fastq = key_list[key_list.index(srr) + 1]
        print("\n", flush=True)
        print("downloading next fastq " + next_fastq+".fastq ...", flush=True)
        process = subprocess.Popen(["./sratoolkit.3.0.2-ubuntu64/bin/fasterq-dump", "-v", "--concatenate-reads", next_fastq], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except IndexError:
        print("No files left to download\n", flush=True)
        pass
    
    print ("downloading current fastq file " + srr + "\n", flush=True)
    subprocess.call(["./sratoolkit.3.0.2-ubuntu64/bin/fasterq-dump", "-v", "--concatenate-reads", srr])

    #rename the current fastq file
    current_name = srr + ".fastq"
    new_name = dictionary[srr]+"_"+current_name
    os.rename(current_name, new_name)
    print ("renaming fastq " + current_name + " to " + new_name + "\n", flush=True)
    
    #find the gfa file to update
    gfa_path = "/home/rmaruzani/storage/vus/toolkits/bifrost"
    gfa = []
    for filename in os.listdir(gfa_path):
        if filename.endswith('.gfa'):
            gfa.append(filename)
            print("Bifrost graph found, called " + gfa[0], flush=True)
        else:
            pass

    #store graph and colours name        
    prefix = os.path.splitext(gfa[0])[0]
    bfg_colors = prefix + ".bfg_colors"
    old_gfa = gfa[0]
    
    #run Bifrost update 
    with open('bifrost_output.log', 'a') as log_file:
        print ("updating Bifrost index with the fastq " + new_name , flush=True)
        subprocess.call(["Bifrost", "update", "-v", "-t", "60", "-k", "31", "-c", "-i", "-d", "-s", new_name, "-g", gfa[0], "-f", bfg_colors, "-o", srr], stdout=log_file, stderr=log_file)
        #remove old fastq, gfa and bfg_colors files
        print ("removing old fastq ..." + new_name, flush=True)
        subprocess.call(["rm", new_name])
        print ("removing old gfa ..." + old_gfa, flush=True)
        subprocess.call(["rm", old_gfa])
        print ("removing old bfg_colors ..." + old_gfa, flush=True)
        subprocess.call(["rm", bfg_colors])
    #wait for the next fastq to finish downloading
    process.wait()
    print ("next fastq downloaded, called " + next_fastq+".fastq" + " moving on to next accession ...", flush=True)

#remove temporary fasterq directories    
tmp_path = "/home/rmaruzani/storage/vus/toolkits/bifrost"
for d in os.listdir(tmp_path):
    if d.startswith('fasterq.tmp'):
        print("removing temporary fasterq directory " + d, flush=True)
        os.rmdir(os.path.join(tmp_path, d))