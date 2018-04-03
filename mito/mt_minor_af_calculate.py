import pandas as pd
import numpy as np 


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
    return sorted_basecounts[1][0]

def reference_allele(row):
    counts = dict({'A': row['A'], 'C': row['C'], 'G': row['G'], 'T': row['T']})
    allele = '' 
    for key, val in counts.items():
        if val == row['ref']:
            allele=key
    return allele


data = pd.read_table('MT_all_positions.txt', sep='\t', names=['type', 'bam', 'pos', 'ref', 'total', 'A', 'C', 'G', 'T'])

data['minor_af'] = data.apply(minor_af, axis=1)
data['minor_allele'] = data.apply(minor_allele, axis=1)
data['reference_allele'] = data.apply(reference_allele, axis=1)

data.to_csv('MT_all_positions_with_af.txt', sep='\t', index=False)


