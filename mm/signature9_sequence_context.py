#!/home/users/cjyoon/anaconda3/bin/python
import sys
import os
import subprocess
import shlex
import re
import argparse
import pysam
import cyvcf2

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='Input vcf with signature 9 only, that will be used to extract +- 10 bp sequence context')
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory')
    parser.add_argument('-r', '--reference', required=False, default='/home/users/cjyoon/reference/GRCh37/human_g1k_v37.fasta')
    args = vars(parser.parse_args())

    return args['input'], args['reference'], args['output_dir']

def main():
    input_vcf, reference, output_dir = argument_parser()
    input_base = os.path.basename(input_vcf)
    output_sequences = os.path.join(output_dir, re.sub(r'.vcf$', '.vcf.21bp.fasta', input_base))
    ref = pysam.FastaFile(reference)

    with open(output_sequences, 'w') as f:
        for variant in cyvcf2.VCF(input_vcf):
            f.write(f'\n>{input_base}_{variant.CHROM}_{variant.POS}_{variant.REF}')
            sequence = ref.fetch(variant.CHROM, variant.POS - 11, variant.POS + 10)
            f.write(f'\n{sequence}')


if __name__=='__main__':
    main()


