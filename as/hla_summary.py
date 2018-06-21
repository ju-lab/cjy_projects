# summarise allele count outputs from polysolver run by ayh. 
# ayh's outputfile contained both SC and SH which are not independent
# Thus recalculate so that only proband is included
# output of this file will be combined with HLA population frequency annotated here: /home/users/team_projects/AnkylosingSpondylitis/analysis/polysolver/AS_winner.hla.comp.txt

with open('AS_winner.hla.txt', 'r') as f:
    for line in f:
        if linecount > 6:
            if not line.startswith(':'):
                if not line.endswith('.txt'):
                    if line.startswith('HLA-A'):
                        _, allele_1, allele_2 = line.strip().split()
                        a_count.append(allele_1)
                        a_count.append(allele_2)
                    elif line.startswith('HLA-B'):
                        _, allele_1, allele_2 = line.strip().split()
                        b_count.append(allele_1)
                        b_count.append(allele_2)
                    elif line.startswith('HLA-C'):
                        _, allele_1, allele_2 = line.strip().split()
                        c_count.append(allele_1)
                        c_count.append(allele_2)
                    else:
                        pass
        linecount += 1

with open('/home/users/cjyoon/Projects/AnkylosingSpondylitis/b27neg_hla_allele_count_15exome_and_proband.txt', 'w') as f:
    f.write('hla_type\tcount\n')
    for k, v in Counter(a_count).items():
        f.write(k + '\t' + str(v) + '\n')
    for k, v in Counter(b_count).items():
        f.write(k + '\t' + str(v) + '\n')
    for k, v in Counter(c_count).items():
        f.write(k + '\t' + str(v) + '\n')



