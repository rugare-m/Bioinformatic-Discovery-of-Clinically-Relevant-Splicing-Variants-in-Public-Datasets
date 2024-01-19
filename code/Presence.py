import csv
import pandas as pd
import pysam
import glob
import subprocess

#extract chromosome and position from the query name
def extract_chromosome_position(key):
    chrom, pos, *_ = key.split('|')
    return chrom, pos

#find matching variant in the vcf file
def find_matching_variant(vcf_file, chrom, pos):
    # Open the VCF file for reading
    vcf_reader = pysam.VariantFile(vcf_file)

    # Iterate through each variant in the VCF file
    for record in vcf_reader:
        if record.chrom == chrom and record.pos == pos:
            return record  # Return the matching variant record

    return None  # If no matching variant is found

#func returns a dictionary with the query name(ie SNP) as the key and samples that have query sequence as the value
def process_tsv(filename):
    result = {}
    
    with open(filename, 'r') as tsv_file:
        reader = csv.reader(tsv_file, delimiter='\t')
        header = next(reader)
        
        for row in reader:
            query_name = row[0]
            values = [header[i] for i in range(1, len(row)) if row[i] == '1']
            
            if values:  # Skip rows with no '1' values
                result[query_name] = values
    
    return result

#query files
reference_qry = './tsv/allsnps_ref_queries.tsv'
alternate_qry = './tsv/allsnps_alt_queries.tsv'

# get dictionaries of queries and samples
print("Processing TSV files...", flush=True)
reference_dict = process_tsv(reference_qry)
alternate_dict = process_tsv(alternate_qry)
print("TSV files processed.", flush=True)


# Remove hg38.fa from the values in the dictionaries - causes problems later
value_to_remove = "/data/sra/hydra2/vus/toolkits/bifrost/hg38.fa"

# Filter the dictionaries to remove hg38.fa
print("Filtering dictionaries...", flush=True)
ref_filtered_dict = {key: [value for value in values if value != value_to_remove] for key, values in reference_dict.items()}
alt_filtered_dict = {key: [value for value in values if value != value_to_remove] for key, values in alternate_dict.items()}
print("Dictionaries filtered.", flush=True)

#Rename the keys in the dictionaries to match the VCF file
print("Rename variants...", flush=True)
new_ref_filtered_dict = {}
for key, value in ref_filtered_dict.items():
    new_key = '|'.join(key.split('|')[:2])
    new_ref_filtered_dict[new_key] = value

new_alt_filtered_dict = {}
for key, value in alt_filtered_dict.items():
    new_key = '|'.join(key.split('|')[:2])
    new_alt_filtered_dict[new_key] = value
print("Variants renamed.", flush=True)


# Remove common values from both dictionaries per matching key, and remove keys that have no values
# We only want to keep the studies that are unique to reference or alternate

print("Removing common values...", flush=True)
keys_to_delete = []
for key in new_alt_filtered_dict:
    if key in new_ref_filtered_dict:
        common_values = set(new_alt_filtered_dict[key]) & set(new_ref_filtered_dict[key])
        new_alt_filtered_dict[key] = list(set(new_alt_filtered_dict[key]) - common_values)
        new_ref_filtered_dict[key] = list(set(new_ref_filtered_dict[key]) - common_values)
    if not new_alt_filtered_dict[key]:
        keys_to_delete.append(key)

for key in keys_to_delete:
    if key in new_alt_filtered_dict:
        del new_alt_filtered_dict[key]
    if key in new_ref_filtered_dict:
        del new_ref_filtered_dict[key]
print("Common values removed.", flush=True)


# Filter the dictionaries based on the number of samples per variant
print("Filtering dictionaries based on count...", flush=True)
n_alt = 5
n_ref = 5
alt_filtered_dict_x = {key: value for key, value in new_alt_filtered_dict.items() if len(value) >= n_ref}
ref_filtered_dict_x = {key: value for key, value in new_ref_filtered_dict.items() if len(value) >= n_alt}
print("Dictionaries filtered based on count.", flush=True)


# Get the common keys between the two dictionaries
print("Finding common keys...", flush=True)
common_keys = set(alt_filtered_dict_x.keys()) & set(ref_filtered_dict_x.keys())
print("Common keys found.", flush=True)

# Create new dictionaries with the common keys/variants and their respective samples
# We want variants appearing in both reference and alternate, but with exclusive samples
alt_filtered_dict_y = {key: value for key, value in alt_filtered_dict_x.items() if key in common_keys}
ref_filtered_dict_y = {key: value for key, value in ref_filtered_dict_x.items() if key in common_keys}

##################################################################################################################################################################################################

# Create a dataframe from the dictionary
# Each study with each variant is a row
def create_dataframe(vcf_file, gtf_file, dictionary):
    keys = dictionary.keys()
    dataframe_lists = []

    # Read GTF file and extract gene start, end positions, and chromosome
    gtf_data = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)
    #gtf_data = gtf_data[gtf_data[2] == 'exon']
    gene_positions = {}
    for _, row in gtf_data.iterrows():
        attributes = row[8]
        gene_id = attributes.split(';')[0].split(' ')[1].strip('"')
        gene_start = int(row[3])
        gene_end = int(row[4])
        chromosome = row[0]
        
        if gene_id in gene_positions:
            if gene_start < gene_positions[gene_id]['start']:
                gene_positions[gene_id]['start'] = gene_start
            if gene_end > gene_positions[gene_id]['end']:
                gene_positions[gene_id]['end'] = gene_end
        else:
            gene_positions[gene_id] = {'start': gene_start, 'end': gene_end, 'chr': chromosome}

    for key in keys:
        variant = []
        chrom, pos = extract_chromosome_position(key)
        variant.append(chrom)
        variant.append(int(pos))

        matching_variant = find_matching_variant(vcf_file, variant[0], variant[1])
        gene_info = matching_variant.info.get('GENEINFO')
        variant_id = matching_variant.id
        if gene_info:
            gene_symbol = gene_info.split('|')[0].split(':')[0]
            if gene_symbol in gene_positions:
                gene_start = gene_positions[gene_symbol]['start']
                gene_end = gene_positions[gene_symbol]['end']
                chromosome = gene_positions[gene_symbol]['chr']
            else:
                gene_start = None
                gene_end = None
                chromosome = None
        else:
            gene_symbol = None
            gene_start = None
            gene_end = None
            chromosome = None

        for val in dictionary[key]:
            row = []
            srr_ids = []
            srr_id = val.split('-')[2].split('.')[0]
            row.append(srr_id)
            row.append(gene_symbol)
            row.append(chromosome)
            row.append(gene_start)
            row.append(gene_end)
            row.append(variant_id)
            dataframe_lists.append(row)

    df = pd.DataFrame(dataframe_lists)
    column_names = ['Study','Gene', 'Chr', 'Gene Start', 'Gene Stop', 'Variant ID']
    df.columns = column_names

    return df


# Define the VCF file
vcf_file = '/home/rmaruzani/storage/hydra2/vus/data/fasta_tools/snps.clinvar.filter.vcf'
gff_file = '/home/rmaruzani/storage/hydra2/vus/data/hg38.ncbiRefSeq.gtf'

print("Creating dataframes...", flush=True)
df_alt = create_dataframe(vcf_file, gff_file, alt_filtered_dict_y)
df_sorted_alt = df_alt.sort_values(by=df_alt.columns[1])
df_ref = create_dataframe(vcf_file, gff_file, ref_filtered_dict_y)
df_sorted_ref = df_ref.sort_values(by=df_ref.columns[1])
print("Dataframes created.", flush=True)

df_sorted_ref.to_csv('reference.csv', index=False)
df_sorted_alt.to_csv('alternate.csv', index=False)

