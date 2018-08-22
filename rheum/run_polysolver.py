'''
This script runs polsolver using python's multiprocessing
Polysolver is installed on the workstation,
This script was run on the workstation with sever mounted
2018.07.05 CJY
'''

import os
import re
import multiprocessing as mp 
import subprocess
import shlex

bamdir = '/home/users/cjyoon/Projects/rheum/bam'
outputdir = '/home/users/cjyoon/Projects/rheum/data_processing/05_polysolver'
bams = [os.path.abspath(os.path.join(bamdir, i)) for i in os.listdir(bamdir) if i.endswith('bam')]
print(bams)

def run_polysolver(bampath):
    samplename = re.sub('.bam$', '', os.path.basename(bampath))
    print(samplename)
    cmd = f'sh /home/users/kjyi/bin/polysolver --bam {bampath} --outdir {outputdir}/{samplename}'
    print(cmd)
    execution = subprocess.Popen(shlex.split(cmd))
    execution.wait()

    return 0 

with mp.Pool(processes=10) as pool:
     pool.map(run_polysolver, bams)



