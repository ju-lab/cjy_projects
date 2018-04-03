import subprocess
import shlex
import os
import sys
import re
local_pcawg_bam_list = '/home/users/cjyoon/Projects/pcawg/md5sum_result.txt'

output_dir = '/home/users/cjyoon/Projects/mito/bam'
with open(local_pcawg_bam_list, 'r') as f:
    for line in f:
        path, status = line.strip().split()
        if os.path.isfile(path) and status == 'True':
            mtbam = os.path.join(output_dir, re.sub(string=os.path.basename(path), pattern=r'.bam$', repl='.MT.bam'))
            print(mtbam)
            cmd = f'samtools view -hbo {mtbam}  {path} MT'
            extract = subprocess.Popen(shlex.split(cmd))
            extract.wait() 
            cmd = f'samtools index {mtbam}'
            index = subprocess.Popen(shlex.split(cmd))
            index.wait() 
        



