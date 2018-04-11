import os
import subprocess
import re
import shlex


filter_directory = '/home/users/cjyoon/Projects/myeloma/analysis/snv_filter'
nin_filter_threshold = 0.01
with open('/home/users/cjyoon/Projects/myeloma/sample_info/tumor_match_abs_all.txt', 'r') as f:
    for line in f:
        sampleID, tumorBam, normalBam = line.strip().split()
        tumorID = re.sub(string=re.sub(string=tumorBam, pattern=r'/home/users/cjyoon/Projects/myeloma/bam/', repl=''), pattern=r'.sorted.md.indel.br.bam', repl='')
        normalID = re.sub(string=re.sub(string=normalBam, pattern=r'/home/users/cjyoon/Projects/myeloma/bam/', repl=''), pattern=r'.sorted.md.indel.br.bam', repl='')
        mutect_vcf = f'/home/users/cjyoon/Projects/myeloma/analysis/mutect/{sampleID}.mutect.vcf'
        strelka_vcf = f'/home/users/cjyoon/Projects/myeloma/analysis/strelka/{sampleID}_{tumorID}_{normalID}/results/variants/somatic.snvs.passonly.vep.vcf.gz'
#        if os.path.isfile(mutect_vcf) and os.path.isfile(manta_vcf) and os.path.isfile(delly_vcf):
        
        # Step 1
        # Apply Not In Normal Filter for Mutect
        cmd = f'python /home/users/cjyoon/scripts/variantFilter/snv/mutect_ninFilter.py -m {mutect_vcf} -s {sampleID} -t {nin_filter_threshold} -o {filter_directory}' 
        print(cmd)
        subprocess.call(shlex.split(cmd))

        filtered_mutect = os.path.join(filter_directory, re.sub(string=os.path.basename(mutect_vcf), pattern=r'.vcf$', repl='.nin_filter_' + str(nin_filter_threshold) + '.vcf'))
        print(filtered_mutect)
        
        # Step 2
        # Intersect filtered Mutect with Strelka
        cmd = f'python /home/users/cjyoon/scripts/variantFilter/snv/intersect_mutect_strelka.py --mutect {filtered_mutect} --strelka {strelka_vcf} -s {sampleID} -o {filter_directory}'
        print(cmd)
        subprocess.call(shlex.split(cmd))
        
        intersect_vcf = os.path.join(filter_directory, sampleID + '_mns.vcf') 
        # Step 3
        # Annotate Filtered VCF with VEP
        cmd = f'python /home/users/cjyoon/scripts/annotate_vcf/annotate_vcf.py -i {intersect_vcf} -o {filter_directory}'
        print(cmd)
        subprocess.call(shlex.split(cmd))
        

