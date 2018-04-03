import subprocess
import shlex


resubmit_list = []
with open('need_rerun.txt', 'r') as f:
	for line in f:
		resubmit_list.append(line.strip())

for i in resubmit_list:
	cmd = "bsub -G team78-grp -P team78-grp -n1 -R'span[ptile=1]' -M24000 -R'select[mem>24000] rusage[mem=24000]' -q long -e ./log3/pos" + str(i) + ".err -o ./log3/pos" + str(i) + ".out python mt_background_vaf.py " + str(i)
	print(cmd)

	execute = subprocess.call(shlex.split(cmd))


