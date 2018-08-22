#!/home/users/cjyoon/anaconda3/bin/python
'''
Script to extract those with T>G mutations with TpTpT, TpTpA, CpTpT context froma given VCF file
2018.07.11 CJY
'''
import sys
import os
import subprocess
import shlex
import re
import argparse
import cyvcf2
import pysam

class Position():
    ''' python class for handling genomic positions
    0-based
    '''
    def __init__(self, chromosome, start, end, is_bp=None, clipped_reads=None):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.is_bp = is_bp
        self.clipped_reads = None

    def __repr__(self):
        return str(self.chromosome + ":" + str(self.start) + '-' + str(self.end))

    def __str__(self):
        return str(self.chromosome + ":" + str(self.start) + '-' + str(self.end))

    def __len__(self):
        return int(self.end - self.start)

    def __iter__(self):
        for i in range(self.start, self.end):
            yield Position(self.chromosome, i, i+1)

    def __next__(self):
        self.start += 1
        self.end += 1
        return self

    def __hash__(self):
        return hash((self.chromosome, self.start, self.end))

    def __eq__(self, other):
        if isinstance(other, Position):
            return self.chromosome == other.chromosome and self.start == other.start and self.end == other.end
        else:
            print("Not of the same class, cannot compare equality")
            return None

    @classmethod
    def fromstring(cls, position_string):
        if isinstance(position_string, str):
            chromosome = position_string.split(':')[0]
            start = int(position_string.split(':')[1].split('-')[0])
            end = int(position_string.split(':')[1].split('-')[1])
            return Position(chromosome, start, end)
        elif isinstance(position_string, Position):
            return position_string
        else:
            print('position_string has to be either a string class or Position class')
            raise TypeError

    @staticmethod
    def check_format(genomicPositionString):
        if re.search(r'[A-Za-z1-9]+:[0-9]+-[0-9]+', genomicPositionString):
            return True
        else:
            return False

    @staticmethod
    def overlap(position1, position2):
        '''true if position1 and position2 has more than 1 overlapping base'''
        try:
            if isinstance(position1, Position) and isinstance(position2, Position):
                if position1.chromosome == position2.chromosome:
                    if min(position1.end, position2.end) > max (position1.start, position2.start):
                        return True
                    else:
                        return False
                else:
                    return False  # cannot compare if two positions are in different chromosome
            else:
                return None # has to be Posiiton class.
        except:
            Exception

    def extend(self, direction, basepairs):
        """extends objects in by specified base pairs, either upstream, downstream, or both"""
        if direction=="up":
            return Position(self.chromosome, max(0, self.start-basepairs), end)
        elif direction=="down":
            return Position(self.chromosome, self.start, self.end + basepairs)
        elif direction=="both":
            return Position(self.chromosome, max(0, self.start - basepairs), self.end + basepairs)
        else:
            print('direction has to be either up, down, or both')
            raise ValueError

def reverse_complement(seq):
    for base in seq:
        if base not in 'ATCGatcg':
            print("Error: NOT a DNA sequence")
            return None
    seq1 = 'ATCGTAGCatcgtagc'
    seq_dict = {seq1[i]: seq1[i+4] for i in range(16) if i < 4 or 8 <= i < 12}
    return "".join([seq_dict[base] for base in reversed(seq)])

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_vcf', required=True, help='Input VCF to extract signature 9 from')
    parser.add_argument('-r', '--reference', default='/home/users/cjyoon/reference/GRCh37/human_g1k_v37.fasta', help='Reference fasta')
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory')
    args = vars(parser.parse_args())

    return args['input_vcf'], args['reference'], args['output_dir']

def get_trinucleotide(query_position, refbase, altbase, reference):
	ref = pysam.FastaFile(reference)
	variant_base = ref.fetch(region=str(query_position))
	trinucleotide = ref.fetch(region=str(query_position.extend('both', 1)))
	if refbase in ['T', 'C']:
		return refbase, altbase, trinucleotide
	else:
		return reverse_complement(refbase), reverse_complement(altbase), reverse_complement(trinucleotide)

def main():
	input_vcf, reference, output_dir = argument_parser()

	output_vcf = os.path.join(output_dir, re.sub('.vcf$', '.sig9.vcf', os.path.basename(input_vcf)))
	vcf_handle = cyvcf2.VCF(input_vcf)
	output_vcf_handle = cyvcf2.Writer(output_vcf, vcf_handle)

	for variant in cyvcf2.VCF(input_vcf):
		var_position = Position(variant.CHROM, variant.POS, variant.POS)
		refbase, altbase, var_trinucleotide = get_trinucleotide(var_position, variant.REF, variant.ALT[0], reference)

		if var_trinucleotide in ['TTT', 'TTA', 'CTT'] and altbase == 'G':
			output_vcf_handle.write_record(variant)




if __name__=='__main__':
    main()


