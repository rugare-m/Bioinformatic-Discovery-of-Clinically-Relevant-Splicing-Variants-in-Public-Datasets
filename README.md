# Bioinformatic Discovery of Clinically Relevant Splicing Variants in Public Datasets

### Background
This work was done as part of a 3-month internship at Kings College London, and [Guy's & St Thomas Hospital's Genomics](https://www.guysandstthomas.nhs.uk) Innovation Unit (GIU) in London. The project was supervised by [Dr Ali Raza Awan](https://www.linkedin.com/in/ali-awan-phd-51041860/?originalSubdomain=uk), the GIU lead.

### Abstract
Approximately 80% of rare diseases have a genetic component. Genetic testing is routinely used to try and diagnose rare diseases. Variants of Uncertain Significance are variants of which there is conflicting, or no evidence for pathogenic or benign classification. These variants make up approximately 50% of all variants described in the ClinVar database, and they can cause significant distress to patients carrying them. Resolving VUSs may hold information that can help diagnose rare diseases and improve patient care. The aim of this project was to develop a computational pipeline to characterise the effects that VUSs have on transcript splicing isoforms and expression levels. We also aimed to benchmark the pipeline against other established tools for predicting splice variants, including SpliceAI. To build our pipeline, we made use of the RNA sequencing data in the Sequence Read Archive (SRA). I built a local searchable database of RNA sequencing data from the SRA using Bifrost. For each synonymous VUS in ClinVar, I queried the database for samples with the VUS, and samples with the reference allele. The aim was to investigate differences in gene expression and splice patterns between the reference and alternate samples. From 2 Tb of RNA FASTQ files, 15 out of 1233 VUSs had enough reference and alternate samples for differential analysis. The project ran out of time before investigating differences in isoforms. Building a local searchable index was expensive in terms of storage and compute time power. 

### Pipeline Workflow 
![VUS](https://github.com/rugare-m/Bioinformatic-Discovery-of-Clinically-Relevant-Splicing-Variants-in-Public-Datasets/assets/88198662/737c29df-f530-44e7-ba59-fcd3f8277f46)

### Collecting RNA Sequencing metadata from Sequence Read Archive
```
BioEntrez.py & RunSelector.py
```

MetaSRA is a webtool that annotates SRA data with tissue of origin information using SRA metadata. I downloaded a list of SRA Study IDs based on 23 tissues adopted from tissues described in the GTEx database. SRA Study IDs hold individual sample Run IDs from a single experiment. SRA currently only supports downloading data using Run IDs, so I had to extract Run ID information for each Study ID, and for each tissue. Additionally, we were only interested in human samples with public access consent. Using MetaSRA Study IDs and the National Center for Biotechnology Information (NCBI)’s Entrez database system, I collected sample Run IDs, Sample Size (Mb, Gb), organism, and consent status. 

### Downloading SRA data and building local searchable RNA Sequencing Database
```
SRAParallelDownload.py & Presence.py & FlankSeq.py
```

I downloaded up to 20 Gb of RNA sequencing samples per tissue to build the local database using the SRA Toolkit. For paired end data, forward and reverse reads were concatenated into a single FASTQ file. The size metadata in SRA the apparently vastly underestimates the actual size of files - so requesting 20 Gb, multiplied by 23 tissues gave approximately 2 Tb of compressed (bgzip) data, not 460 Gb. Next, I indexed the approximately 2 Tb of FASTQ files in the local database enabling it to be searchable. For querying the database, I used Bifrost - a tool that indexes sequence data by building compacted, coloured de Buijn graphs. The Bifrost index is searchable by inputting a FASTA file, with the output being a CSV file with database samples containing each sequence in the FASTA file. To generate the alternate (VUS) query file, I created a FASTA file where each sequence was a VUS nucleotide flanked by 15 reference nucleotides. A reference query file was generated with the same method, except the VUS was replaced with the reference allele.  

### Quality control and differential expression analysis
```
GeneReadsExport.nf
```

VUSs that had at least 5 samples supporting the reference and alternate alleles were carried forward to investigate the differences in expression and splicing. Before this, some sanity checks were performed on data. First, I checked to see if the query sequence in each positive sample was found in the expected region. Bifrost will return a match for a query sequence that occurs anywhere within the sample – we were only interested in matches at the loci of the VUS. I mapped the FASTQ files to GRCh38, using the STAR aligner and exported reads mapping to the VUS gene. Next, I used the grep tool in bash to confirm the query sequence was still present in the gene only reads – if, indeed there were still reads present. Out of 1233 synonymous VUS SNPs, 15 passed the quality control and sanity check. I first performed differential expression analysis on these samples using Salmon, which showed no differences in gene expression. At this point we ran out of time on the project before investigating any differences in transcript splicing between reference and VUS allele samples. 

### Thanks to
