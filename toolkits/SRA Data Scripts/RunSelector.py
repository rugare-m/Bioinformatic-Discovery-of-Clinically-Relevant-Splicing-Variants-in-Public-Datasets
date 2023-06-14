import csv
import pandas as pd
import random
import argparse

'''
usage: python RunSelector.py -s 51200 -c tissue_batch.csv -i recount3_x_metasra.txt -ot download_runs.txt -oc download_runs.csv
-c is the output csv file from BioEntrez.py
-i is a txt of Study IDs, can be obtained from column one of output csv file from BioEntrez.py
-ot is the output txt file, contains a list of study IDs
-oc is the output csv file, contains the dataframe of study id, run id, consent, size, organism, tissue, for each run
'''
parser = argparse.ArgumentParser(
                    prog='run_selector.py',
                    description='select runs from a csv file based on a study id list, there is a size limit for runs in each study')

parser.add_argument('-c', '--csv')
parser.add_argument('-ot', '--out_txt')
parser.add_argument('-oc', '--out_csv')
parser.add_argument('-i', '--sid', type=str, help='study id')
parser.add_argument('-s', '--size', type=int)
args = parser.parse_args()

#open the csv file and convert to dictionary
df = pd.read_csv(args.csv, index_col=False)
sort_df = df.sort_values(['Study', 'Size (MB)'])

download_list = []
input_txt = args.sid

with open(input_txt, 'r') as f:
    my_list = f.read().splitlines()

#select runs from the csv file based on the study id list
for item in my_list:
    size = 0
    for index, row in sort_df.iterrows():
        if row['Study'] == item:
            accession = []
            if size >= args.size:
                print ("Size limit exceeded for study" +" "+ str(item) +", "+ "removing last accession "+str(row['Run ID']) + " and proceeding to next study ...")
                del download_list[-1]
                break
            accession.append(row['Run ID'])
            accession.append(row['Size (MB)'])
            accession.append(row['Tissue'])
            accession.append(row['Study'])
            #assigning 5 letter tissue codes
            if row['Tissue'] == 'brain':
                accession.append("BRAIN")
            elif row['Tissue'] == 'adipose':
                accession.append("ADPSE")
            elif row['Tissue'] == 'adrenal gland':
                accession.append("ADGLD")
            elif row['Tissue'] == 'aorta':
                accession.append("AORTA")
            elif row['Tissue'] == 'blood plasma':
                accession.append("BPLSM")
            elif row['Tissue'] == 'breast':
                accession.append("BREST")
            elif row['Tissue'] == 'cardiac atrium':
                accession.append("CARAT")
            elif row['Tissue'] == 'colon':
                accession.append("COLON")
            elif row['Tissue'] == 'left ventricle':
                accession.append("LVENT")
            elif row['Tissue'] == 'liver':
                accession.append("LIVER")
            elif row['Tissue'] == 'lung':
                accession.append("LUNG")
            elif row['Tissue'] == 'lymphocyte':
                accession.append("LYMPH")
            elif row['Tissue'] == 'Oesophagus mucosa':
                accession.append("MUCOS")
            elif row['Tissue'] == 'ovary':
                accession.append("OVARY")
            elif row['Tissue'] == 'pancreas':
                accession.append("PANCR")
            elif row['Tissue'] == 'prostate gland':
                accession.append("PROGL")
            elif row['Tissue'] == 'skeletal muscle':
                accession.append("SKMUS")
            elif row['Tissue'] == 'small intestine':
                accession.append("INTES")
            elif row['Tissue'] == 'spleen':
                accession.append("SPLEE")
            elif row['Tissue'] == 'stomach':
                accession.append("STOMA")
            elif row['Tissue'] == 'testis':
                accession.append("TESTI")
            elif row['Tissue'] == 'tibia':
                accession.append("TIBIA")
            elif row['Tissue'] == 'uterus':
                accession.append("UTERI")
            
                
            download_list.append(accession)
            size += row['Size (MB)']

download_df = pd.DataFrame(download_list, columns =["Run ID","Size (MB)", "Tissue", "Study", "Tissue Code"])
run_ids = download_df['Run ID']  

text_outfile = args.out_txt
csv_outfile = args.out_csv

print ("Generating output files...")
download_df.to_csv(csv_outfile, index=False)
run_ids.to_csv(text_outfile, index=False, header=False, line_terminator='\n')


