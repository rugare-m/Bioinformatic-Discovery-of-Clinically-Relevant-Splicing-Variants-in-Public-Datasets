import pysam
from Bio import SeqIO
from fasta_reader import read_fasta
from fasta import FASTA
from cyvcf2 import VCF, Writer
import pyfastx
import argparse
import random

parser = argparse.ArgumentParser(
                    prog='FlankSeq',
                    description='FlankSeq will take in a VCF and reference genome and output a FASTA of variants described in the VCF flanked by sequences of a user defined length')

parser.add_argument('-v', '--vcf')
parser.add_argument('-r', '--reference')
parser.add_argument('-o', '--output')
parser.add_argument('-f', '--flank_length', type=int)
parser.add_argument('-s', '--snp')

args = parser.parse_args()


hg38 = pyfastx.Fasta(args.reference)


snp_filter_clinvar = VCF(args.vcf)

f = open(args.output, "x")

snp = args.snp.upper()

for variant in snp_filter_clinvar:
    if snp == "REF":
        if variant.CHROM[3].isdigit():
            flank_length = (int(args.flank_length)-1)
            intervalp = (variant.POS, variant.POS+flank_length)
            intervaln = (variant.POS-flank_length, variant.POS)
            fasta = []
            fasta.append(hg38.fetch(variant.CHROM,intervalp))
            fasta.append(variant.REF)
            fasta.append(hg38.fetch(variant.CHROM,intervaln))
            seq = []
            seq.append("".join(fasta))
            head = []
            head.append(">")
            head.append(variant.CHROM)
            head.append("|")
            head.append(str(variant.POS))
            head.append("|")
            head.append("REF")
            head.append("|")
            head.append(variant.REF)
            header = []
            header.append("".join(head))
            #print("writing read ...")
            f.write(str(header[0])+"\n"+seq[0]+"\n")
    
    elif snp == "ALT":         
        if variant.CHROM[3].isdigit():
            flank_length = (int(args.flank_length)-1)
            
            intervalp = (variant.POS, variant.POS+flank_length)
            intervaln = (variant.POS-flank_length, variant.POS)
            fasta = []
            fasta.append(hg38.fetch(variant.CHROM,intervalp))
            fasta.append(variant.ALT[0])
            fasta.append(hg38.fetch(variant.CHROM,intervaln))
            seq = []
            #join the fasta seq
            seq.append("".join(map(str, fasta)))
            #print (seq)
            head = []
            head.append(">")
            head.append(variant.CHROM)
            head.append("|")
            head.append(str(variant.POS))
            head.append("|")
            head.append("ALT")
            head.append("|")
            head.append(variant.ALT[0])
            header = []
            header.append("".join(map(str, head)))
            #print("writing sequence ...")
            f.write(str(header[0])+"\n"+seq[0]+"\n")
f.close()


