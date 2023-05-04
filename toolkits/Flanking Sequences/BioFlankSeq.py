 #!/usr/bin/env python3
 # -*- coding: utf-8 -*-

import pysam
from Bio import SeqIO
import numpy
from fasta_reader import read_fasta
from fasta import FASTA
from cyvcf2 import VCF, Writer
import pyfastx
import argparse
import random

def get_flanking_sequences(fasta_file, sequence_id, position, flank_length):
    """
    Args:
        fasta_file (str): Path to the FASTA file.
        sequence_id (str): ID of the sequence in the FASTA file.
        position (int): Position of the desired flanking sequences.
        flank_length (int): Length of the flanking sequences.
    """
    left_flank = ""
    right_flank = ""

    # Parse the FASTA file and retrieve the desired sequence
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == sequence_id:
            sequence = record.seq
            if position >= flank_length and position + flank_length < len(sequence):
                left_flank = sequence[position - flank_length:position]
                right_flank = sequence[position:position + flank_length]
                break

    return left_flank, right_flank


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
            flank_length = (int(args.flank_length))
            intervalp = (variant.POS, variant.POS+flank_length)
            intervaln = (variant.POS-flank_length, variant.POS)
            fasta = []
            
            left_flank, right_flank = get_flanking_sequences(args.reference, variant.CHROM, variant.POS, flank_length)
            
            
            fasta.append(str(left_flank.upper()))
            #print (left_flank)
            fasta.append(variant.REF)
            fasta.append(str(right_flank.upper()))
            
        
            seq = []
            seq.append("".join(fasta))
            
            #print (seq)
            
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
            print("writing records ...")
            #print(seq)
            f.write(str(header[0])+"\n"+seq[0]+"\n")
    
    elif snp == "ALT":         
        if variant.CHROM[3].isdigit():
            flank_length = (int(args.flank_length))
            
            intervalp = (variant.POS, variant.POS+flank_length)
            intervaln = (variant.POS-flank_length, variant.POS)
            fasta = []
            
            left_flank, right_flank = get_flanking_sequences(args.reference, variant.CHROM, variant.POS, flank_length)
            
            fasta.append(str(left_flank.upper()))
            #print (left_flank)
            fasta.append(variant.ALT[0])
            fasta.append(str(right_flank.upper()))
            
            
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
            print("writing sequence ...")
            f.write(str(header[0])+"\n"+seq[0]+"\n")
f.close()
