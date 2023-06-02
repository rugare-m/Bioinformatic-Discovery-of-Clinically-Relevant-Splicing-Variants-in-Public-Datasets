#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
import pyfastx
import argparse
from cyvcf2 import VCF


def get_flanking_sequences(fasta, sequence_id, position, flank_length):
    """
    Args:
        fasta (pyfastx.Fasta): Fasta object containing the reference genome.
        sequence_id (str): ID of the sequence in the FASTA file.
        position (int): Position of the desired flanking sequences.
        flank_length (int): Length of the flanking sequences.
    """
    left_flank = ""
    right_flank = ""

    sequence = fasta[sequence_id]

    if position >= flank_length and position + flank_length < len(sequence):
        lpos = position - 1
        rpos = position
        left_flank = sequence[lpos - flank_length:lpos]
        right_flank = sequence[rpos:rpos + flank_length]

    return left_flank, right_flank


parser = argparse.ArgumentParser(
    prog='FlankSeq',
    description='FlankSeq will take in a VCF and reference genome and output a FASTA of variants described in the VCF flanked by sequences of a user defined length'
)

parser.add_argument('-v', '--vcf')
parser.add_argument('-r', '--reference')
parser.add_argument('-o', '--output')
parser.add_argument('-f', '--flank_length', type=int)
parser.add_argument('-s', '--snp')

args = parser.parse_args()

hg38 = pyfastx.Fasta(args.reference)

results = []

with open(args.output, "w") as f:
    snp = args.snp.upper()

    for variant in VCF(args.vcf):
        if variant.CHROM:
            flank_length = int(args.flank_length)
            gapp = variant.POS + 1
            gapn = variant.POS - 1
            intervalp = (gapp, variant.POS + flank_length)
            intervaln = (variant.POS - flank_length, variant.POS)
            fasta = []

            left_flank, right_flank = get_flanking_sequences(
                hg38, variant.CHROM, variant.POS, flank_length
            )

            fasta.append(str(left_flank))
            fasta.append(str(variant.REF if snp == "REF" else variant.ALT[0]))
            fasta.append(str(right_flank))

            seq = "".join(fasta)

            head = ">{}|{}|{}|{}".format(
                variant.CHROM,
                variant.POS,
                "REF" if snp == "REF" else "ALT",
                variant.REF if snp == "REF" else variant.ALT[0]
            )

            results.append("{}\n{}\n".format(head, seq))

    f.writelines(results)

print("Writing records complete.")
