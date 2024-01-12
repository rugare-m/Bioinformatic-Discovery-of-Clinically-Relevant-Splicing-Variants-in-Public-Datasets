# Bioinformatic Discovery of Clinically Relevant Splicing Variants in Public Datasets

### Background
This work was done as part of a 3-month internship at Kings College London, and [Guy's & St Thomas Hospital's Genomics](https://www.guysandstthomas.nhs.uk) Innovation Unit (GIU) in London. The project was supervised by [Dr Ali Raza Awan](https://www.linkedin.com/in/ali-awan-phd-51041860/?originalSubdomain=uk), the GIU lead.


### Abstract
Approximately 80% of rare diseases have a genetic component. Genetic testing is routinely used to try and diagnose rare diseases. Variants of Uncertain Significance are variants of which there is conflicting, or no evidence for pathogenic or benign classification. These variants make up approximately 50% of all variants described in the ClinVar database, and they can cause significant distress to patients carrying them. Resolving VUSs may hold information that can help diagnose rare diseases and improve patient care. The aim of this project was to develop a computational pipeline to characterise the effects that VUSs have on transcript splicing isoforms and expression levels. To do this we made use of the RNA sequencing data in the Sequence Read Archive (SRA). I built a local searchable database of RNA sequencing data from the SRA using Bifrost. For each synonymous VUS in ClinVar, I queried the database for samples with the VUS, and samples with the reference allele. The aim was to investigate differences in gene expression and splice patterns between the reference and alternate samples. From 2 Tb of RNA FASTQ files, 15 out of 1233 VUSs had enough reference and alternate samples for differential analysis. Building a local searchable index was expensive in terms of storage and compute time power. 

### Pipeline Workflow 

### Bulk downloading RNA Seq data from the Sequence Read Archive 

### Building local searchable RNA Seq Database

### Building Query FASTA files

### Preparing reads for expression and alternative splicing analysis
