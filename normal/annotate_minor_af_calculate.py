import pandas as pd
import numpy as np 
import os
import argparse


def minor_af(row):
    counts = dict({'A': row['A'], 'C': row['C'], 'G': row['G'], 'T': row['T']})
    sorted_basecounts = ([(k, v) for k,v in sorted(counts.items(), key=lambda item: item[1], reverse=True) ])
    if int(row['total']) != 0:
        return int(sorted_basecounts[1][1])/int(row['total'])
    else:
        return np.nan

def minor_allele(row):
    counts = dict({'A': row['A'], 'C': row['C'], 'G': row['G'], 'T': row['T']})
    sorted_basecounts = ([(k, v) for k,v in sorted(counts.items(), key=lambda item: item[1], reverse=True) ])
    if sorted_basecounts[1][1]==0 and sorted_basecounts[2][1]==0 and sorted_basecounts[3][1]==0:
        # if all bases except the major base count is 0 then there is no minor allele 
        return 'N'
    else:
        return sorted_basecounts[1][0]

def reference_allele(row):
    counts = dict({'A': row['A'], 'C': row['C'], 'G': row['G'], 'T': row['T']})
    allele = '' 
    for key, val in counts.items():
        if val == row['ref']:
            allele=key
    return allele

def argument_parser():
    parser = argparse.ArgumentParser(description='Add minor allele fraction, REF, ALT base as columns to the count table')
    parser.add_argument('-i', '--input', required=True, help='Input count table')
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory')

    args = vars(parser.parse_args())
    
    return args['input'], args['output_dir']

def main():
    input, output_dir = argument_parser()

    data = pd.read_table(input, sep='\t', names=['bam', 'chromosome', 'position', 'ref', 'total', 'A', 'C', 'G', 'T', 'position_name'])

    data['minor_af'] = data.apply(minor_af, axis=1)
    data['minor_allele'] = data.apply(minor_allele, axis=1)
    data['reference_allele'] = data.apply(reference_allele, axis=1)
    outputpath = os.path.join(output_dir, os.path.basename(input) + '_with_af.txt')

    data.to_csv(outputpath, sep='\t', index=False)

if __name__ == '__main__':
    main()


