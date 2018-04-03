import pysam
import re
mt = pysam.FastaFile('./rCRS/chrMT.fa')
acccc_pos = []
for i in re.finditer(string=mt.fetch('MT'), pattern=r'ACCCC'):
    acccc_pos.append(i.start() + 1)

tcccc_pos = []
for j in re.finditer(string=mt.fetch('MT'), pattern=r'TCCCC'):
    tcccc_pos.append(j.start()+1)


# 10x greater vat rate in these positions=
high_vaf_pos = [72,152,204,248,263,302,310,316,515,567,750,1438,2487,3109,3492,6419,7028,10277,10306,12684,12705,12825,13062,13095,13105,13650,15326,16129,16183,16189,16218,16230,16249,16259,16263,16264,16274,16278,16284,16288,16293,16301,16311,16355,16356,16368,16390,16399,16427,16444,16496,16519,16527]

set(high_vaf_pos) & set(tcccc_pos)
set(high_vaf_pos) & set(acccc_pos)
"""
In [29]: len(acccc_pos)
Out[29]: 79
In [30]: len(tcccc_pos)
Out[30]: 50
In [31]: len(high_vaf_pos)
Out[31]: 53
In [42]: set(high_vaf_pos) & set(acccc_pos)
Out[42]: {302, 567, 2487, 6419, 10277, 16183}
In [43]: set(high_vaf_pos) & set(tcccc_pos)
Out[43]: {310, 16189}
"""
