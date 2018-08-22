"""Submits Sequenza Prep to qsub
March 16 2018 Jongsoo Yoon"""

import subprocess, shlex
reference = '/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta'
with open('/home/users/cjyoon/Projects/myeloma/sample_info/tumor_match_abs_batch2_20180706.txt', 'r') as f:
    for line in f:
        if not line.startswith('#'):
            sampleName, tumorBam, normalBam = line.strip().split()
            cmd = f'qsubmitter.py -n bnode3 -i {sampleName}.seqz.sh -p 2 -c "/home/users/cjyoon/scripts/run_sequenza/sequenza_fromBam.sh {normalBam} {tumorBam} {sampleName} {reference}"'
            print(cmd)
            subprocess.call(shlex.split(cmd))

