import pandas as pd
import re
from Bio import Entrez
Entrez.email = "rugarem@liverpool.ac.uk"

lists = []

runcounter = 0

#doesnt change, this is tissue information from metasra
pd_tissues = pd.read_csv('/home/rmaruzani/storage/vus/toolkits/sra/merged_metasra.csv')
dictionary = dict(zip(pd_tissues['study_id'], pd_tissues['tissue search']))

#take recount3-metasra Study ID instersect, and output "Study", "Run ID", "Consent", "Size (MB)", "Organism", "Tissue" csv and txt files
#human only, no fibroblast
with open('/home/rmaruzani/storage/vus/toolkits/sra/recount3_studies.txt') as f:
    for line in f:
        
        if runcounter >= 100000:
            break
        Entrez.sleep_between_tries = 35
        
        handle = Entrez.esearch(db="sra", term=line, retmax=10000)
        record = Entrez.read(handle)
        run_ids = record['IdList']
        
        #for run_id in run_ids:
        handle = Entrez.efetch(db="sra", id=run_ids, rettype="runinfo")
        run_info = handle.read().decode('utf-8')
                    
        handle.close()
        run_info = run_info.splitlines()
        run_info_header = run_info[0].split(',')
        run_info_data = [dict(zip(run_info_header, line.split(','))) for line in run_info[1:]]
        string = str(run_info_data)

        consent = re.findall(r'<Consent\>(public)\<\/Consent\>', string)
        study = re.findall(r'<SRAStudy\>(\w+)\<\/SRAStudy\>', string)
        size = re.findall(r'<size_MB\>(\d+)\<\/size_MB\>', string)
        run = re.findall(r'<Run\>(\w+)\<\/Run\>', string)
        organism = re.findall(r'<ScientificName\>([A-Z][a-z]+\s[A-Za-z]+)\<\/ScientificName\>', string)
        
        
        zipper = zip(study, run, consent, size,organism)
        mapped = [[a, b, c, d, e] for a, b, c, d, e in zipper]
        
        #append tissue information        
        for item in mapped:
            for key in item:
                if key in dictionary:
                    item.append(dictionary[key])
                    #print(item)

        
        for map in mapped:
            if runcounter >= 100000:
                break
            lists.append(map)
            runcounter += 1
            print (runcounter)
            print(map)
        

#print (lists)
        
df = pd.DataFrame(lists, columns =["Study", "Run ID", "Consent", "Size (MB)", "Organism", "Tissue"])
nofb_df = df[df['Tissue']!="fibroblast"]
nofb_df["Tissue"] = nofb_df["Tissue"].replace(to_replace="subcutaneous adipose tissue", value="adipose")
homo_df = nofb_df[nofb_df['Organism']=="Homo sapiens"]
homo_df.to_csv('tissue_recount3_x_metasra.csv', encoding='utf-8', index=False) 