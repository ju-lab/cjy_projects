#!/home/users/cjyoon/anaconda3/bin/python
import sys
import os
import subprocess
import shlex
import re
import argparse
from collections import Counter

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='')
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory')
    args = vars(parser.parse_args())

    return args['input'], args['output_dir']

def convert_to_hla_nomenclature(polysolver_hla):
    hla_split = polysolver_hla.upper().split('_') 
    rejoin = hla_split[0] + '-' + hla_split[1] + '*' + hla_split[2] + ':' + hla_split[3]
    return rejoin

def main():
    polysolver_dir = '/home/users/cjyoon/Projects/myeloma/analysis/polysolver/'
    sampleNames = [i for i in os.listdir(polysolver_dir) if os.path.isdir(os.path.join(polysolver_dir, i))]
    hla_results = []
    homozygous_present = []
    hla_a = [] 
    hla_b = []
    hla_c = []
    for sample in sampleNames:
        polysolver_output = os.path.join(polysolver_dir, sample + '/' + 'winners.hla.txt')   
        # sample_result = dict() 
        with open(polysolver_output, 'r') as f:
            for line in f:
                hlatype, allele1, allele2 = line.strip().split()
                allele1_nomenclature = convert_to_hla_nomenclature(allele1)
                allele2_nomenclature = convert_to_hla_nomenclature(allele2)
                hla_results.append(allele1_nomenclature)
                hla_results.append(allele2_nomenclature)
                if allele1 == allele2:
                    print(allele1, allele2)
                    homozygous_present.append(sample)

                if hlatype == 'HLA-A':
                    hla_a.append(allele1_nomenclature)
                    hla_a.append(allele2_nomenclature)
                elif hlatype == 'HLA-B':
                    hla_b.append(allele1_nomenclature)
                    hla_b.append(allele2_nomenclature)
                elif hlatype == 'HLA-C':
                    hla_c.append(allele1_nomenclature)
                    hla_c.append(allele2_nomenclature)
                else:
                    raise ValueError

    print(Counter(hla_results))
    print(homozygous_present)
    print(Counter(hla_a))
    print(Counter(hla_b))
    print(Counter(hla_c))

    with open('/home/users/cjyoon/Projects/myeloma/analysis/polysolver/analysis/hla_a_observed_mm.txt', 'w') as f:
        f.write('HLA_allele\tobserved')
        for allele, count in Counter(hla_a).items():
            f.write(f'\n{allele}\t{count}')

    with open('/home/users/cjyoon/Projects/myeloma/analysis/polysolver/analysis/hla_b_observed_mm.txt', 'w') as f:
        f.write('HLA_allele\tobserved')
        for allele, count in Counter(hla_b).items():
            f.write(f'\n{allele}\t{count}')

    with open('/home/users/cjyoon/Projects/myeloma/analysis/polysolver/analysis/hla_c_observed_mm.txt', 'w') as f:
        f.write('HLA_allele\tobserved')
        for allele, count in Counter(hla_c).items():
            f.write(f'\n{allele}\t{count}')


if __name__=='__main__':
    main()


