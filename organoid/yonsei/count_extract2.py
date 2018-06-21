#this script extracts ensembl gene ID, gene length, and count from featureCounts output and converts Ensembl geneID to EntrezID usi ng myGene package. Discards any Ensembl entry without an EntrezID associated. 
# changed to just extracting ensembl ID since mygene querying took so long
# instead using convert2entrez.R for faster conversion using biomart and tidyverse
# 2018.06.13 Jongsoo Yoon

import sys
import re
#import mygene

inputcount = sys.argv[1]
output = re.sub(r'.txt$', '.simple.txt', inputcount)
#mg = mygene.MyGeneInfo()
with open(output, 'w') as g:
    g.write('EnsemblGeneID\tGeneLength\tCount\n')
    with open(inputcount, 'r') as f:
        for line in f:
            if not line.startswith('#') and not line.startswith('Geneid'):
                Geneid, Chr, Start, End, Strand, Length, count = line.strip().split('\t')
                g.write(f'{Geneid}\t{Length}\t{count}\n')
'''
                annot = mg.getgene(Geneid)
                if annot != None:
                    if type(annot) is list:
                        annot = annot[0]
                    else:
                        pass

                    if 'entrezgene' in annot.keys():
                        entrez = annot['entrezgene']
                        g.write(f'{entrez}\t{Length}\t{count}\n')
'''
