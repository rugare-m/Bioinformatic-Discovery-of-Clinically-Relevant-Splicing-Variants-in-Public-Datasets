import csv
import pandas as pd
import random
import argparse

#argparse
parser = argparse.ArgumentParser(
                    prog='run_selector.py',
                    description='get run ids for given tissue')

parser.add_argument('-t', '--tissue')
parser.add_argument('-c', '--csv')
parser.add_argument('-s', '--size', type=int)

args = parser.parse_args()

#open the csv file and convert to dictionary
df = pd.read_csv(args.csv, index_col=False)
dictionary = df.set_index('Run ID').T.to_dict('list')

#func to randomly pick a run id from the dictionary
def pick():       
    return random.choice(list(dictionary.keys()))

#open the tissue txt file
tfile = open(args.tissue,'r')

#and loop through the tissues
for line in tfile:
    line = line.rstrip("\n")
    print (line)
    download_list = []

#group the dataframe by tissue and sum the size of each tissue
    group = df.groupby('Tissue')['Size (MB)'].sum()
    tissue_size = group[line]

    '''
    while the size of the tissue is less than the size specified by the user, 
    randomly pick a run id and add it to the download list, and remove it from the dictionary so it won't be picked again
    if size becomes greater than the size specified by the user, remove the last run id
    and subtract its size from the total size
    '''
    size = 0
    while size <= args.size:
        #in case the size of the tissue is less than the size specified by the user
        if size == tissue_size:
            break
        run = pick()
        mb = int(dictionary[run][2])
        tissue = str(dictionary[run][-1])
        if tissue == line:
            #print (run)
            download_list.append(run)
            size += mb
            del dictionary[run]
            if size > args.size:
                del download_list[-1]
                size -= mb
                break
    print (size)

    #write the run ids to a text file
    output = line + '.txt'
    with open(output, 'w') as fp:
        for srr in download_list:
            fp.write("%s\n" % srr)
        print('Done\n')