import csv
import pandas as pd
import random
import argparse

#argparse
parser = argparse.ArgumentParser(
                    prog='run_selector.py',
                    description='xxx')

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
    # read the contents of the file and split it into a list
    my_list = f.read().splitlines()

print (my_list)

for item in my_list:
    size = 0
    for index, row in sort_df.iterrows():
        if row['Study'] == item:
            accessions = []
            if size >= args.size:
                del download_list[-1]
                break
            accessions.append(row['Run ID'])
            accessions.append(row['Size (MB)'])
            accessions.append(row['Tissue'])
            accessions.append(row['Study'])
            download_list.append(accessions)
            size += row['Size (MB)']

download_df = pd.DataFrame(download_list, columns =["Run ID","Size (MB)", "Tissue", "Study"])
run_ids = download_df['Run ID']  

text_outfile = args.out_txt
csv_outfile = args.out_csv

download_df.to_csv(csv_outfile, index=False)
run_ids.to_csv(text_outfile, index=False, header=False, line_terminator='\n')


