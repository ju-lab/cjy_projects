"""
Dec 27 2017 Jongsoo Yoon (cjyoon@kaist.ac.kr)

This script is used to parallelize delly submission for myeloma cancer/normal pairs. 
Each job is expected to take > 1week

Same job for batch 2 2018.07.06
"""

import subprocess
import re
import sys
import os
import shlex

#sampleFile = '/home/users/cjyoon/Projects/myeloma/sample_info/tumor_match_abs_analysisready.txt'
#sampleFile = '/home/users/cjyoon/Projects/myeloma/sample_info/tumor_match_abs_analysisready_20180123.txt'
sampleFile = '/home/users/cjyoon/Projects/myeloma/sample_info/tumor_match_abs_batch2_20180706.txt' 

with open(sampleFile, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            print(line)
            sampleid, cancerbam, normalbam = line.strip().split()
            cmd =f'qsubmitter.py -i {sampleid}.delly.sh -m 16gb -c "python /home/users/cjyoon/scripts/autobahn/autobahn.py -o /home/users/cjyoon/Projects/myeloma/analysis -i {sampleid} -m {cancerbam} -n {normalbam} -r /home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta -d 0 -l delly"'
            execute = subprocess.Popen(shlex.split(cmd))
            execute.wait()
            print(cmd)


