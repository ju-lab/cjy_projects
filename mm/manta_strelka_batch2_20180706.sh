#!/bin/bash
#PBS -l nodes=bnode1:ppn=4,mem=36gb
#PBS -M julab.job@gmail.com
#PBS -m abe
#PBS -j oe
#PBS -q long
#PBS -o /dev/null
#PBS -e /dev/null
cd $PBS_O_WORKDIR
python /home/users/cjyoon/scripts/autobahn/autobahn.py -o /home/users/cjyoon/Projects/myeloma/analysis -s /home/users/cjyoon/Projects/myeloma/sample_info/tumor_match_abs_batch2_20180706.txt -r /home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta -d 0 -l strelka manta 1> /home/users/cjyoon/Projects/myeloma/analysis/log/manta_strelka_batch2_20180706.sh.stdout 2> /home/users/cjyoon/Projects/myeloma/analysis/log/manta_strelka_batch2_20180706.sh.stderr