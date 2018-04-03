import pysam

mt = pysam.FastaFile('./rCRS/chrMT.fa')

high_background_pos = [72,152,204,248,263,302,310,316,515,567,750,1438,2487,3109,3492,6419,7028,10277,10306,12684,12705,12825,13062,13095,13105,13650,15326,16129,16183,16189,16218,16230,16249,16259,16263,16264,16274,16278,16284,16288,16293,16301,16311,16355,16356,16368,16390,16399,16427,16444,16496,16519,16527]
#for pos in high_background_pos:
#    print(f'>{pos}')
#    print(mt.fetch('MT', pos -5, pos + 5))

with open('high_vaf_flanking_10bp_bothdirection.fa', 'w') as f:
    for pos in high_background_pos:
        f.write(f'>{pos}\n')
        f.write(mt.fetch('MT', pos -5, pos + 5) + '\n')
