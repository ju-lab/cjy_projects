"""Split Tumor and Normal WGS Bam into Separate Rows for easier processing"""

master_table = 'SampleList_101016.txt'
linecount=0
with open('pcawg_korea_data_template.txt', 'w') as g:
    with open(master_table, 'r') as f:
        for line in f:
            columnNumbers = len(line.strip().split('\t'))
            if linecount == 0:
                g.write('final_index\tdonor_unique_id\tdcc_project_code\tsubmitter_donor_id\ticgc_donor_id\twgs_submitter_specimen_id\twgs_icgc_specimen_id\twgs_submitter_sample_id\twgs_icgc_sample_id\twgs_aliquot_id\twgs_alignment_gnos_repo\twgs_alignment_gnos_id\twgs_alignment_bam_file_name\tspecimen\ttumor_specimen_count')
            elif columnNumbers >= 35:
                final_index, our_histology, our_subgroup, donor_unique_id, dcc_project_code, submitter_donor_id, icgc_donor_id, oct2015_donor, santa_cruz_pilot, validation_by_deep_seq, normal_wgs_submitter_specimen_id, normal_wgs_icgc_specimen_id, normal_wgs_submitter_sample_id, normal_wgs_icgc_sample_id, normal_wgs_aliquot_id, normal_wgs_alignment_gnos_repo, normal_wgs_alignment_gnos_id, is_oct2015_normal_wgs_alignment, normal_wgs_alignment_bam_file_name, tumor_wgs_specimen_count, tumor_wgs_submitter_specimen_id, tumor_wgs_icgc_specimen_id, tumor_wgs_submitter_sample_id, tumor_wgs_icgc_sample_id, tumor_wgs_aliquot_id, tumor_wgs_oxog_score, tumor_wgs_alignment_gnos_repo, tumor_wgs_alignment_gnos_id, is_oct2015_tumor_wgs_alignment, tumor_wgs_alignment_bam_file_name, *args = line.strip().split('\t')
                g.write(f'\n{final_index}\t{donor_unique_id}\t{dcc_project_code}\t{submitter_donor_id}\t{icgc_donor_id}\t{normal_wgs_submitter_specimen_id}\t{normal_wgs_icgc_specimen_id}\t{normal_wgs_submitter_sample_id}\t{normal_wgs_icgc_sample_id}\t{normal_wgs_aliquot_id}\t{normal_wgs_alignment_gnos_repo}\t{normal_wgs_alignment_gnos_id}\t{normal_wgs_alignment_bam_file_name}\tnormal\t{tumor_wgs_specimen_count}')
                for i in range(0, int(tumor_wgs_specimen_count)):
                    print(i)
                    tumor_wgs_submitter_specimen_id_i = tumor_wgs_submitter_specimen_id.split(',')[i]
                    tumor_wgs_icgc_specimen_id_i = tumor_wgs_icgc_specimen_id.split(',')[i]
                    tumor_wgs_submitter_sample_id_i = tumor_wgs_submitter_sample_id.split(',')[i]
                    tumor_wgs_icgc_sample_id_i = tumor_wgs_icgc_sample_id.split(',')[i]
                    tumor_wgs_aliquot_id_i = tumor_wgs_aliquot_id.split(',')[i]
                    tumor_wgs_alignment_gnos_repo_i = tumor_wgs_alignment_gnos_repo.split('|')[i]
                    tumor_wgs_alignment_gnos_id_i = tumor_wgs_alignment_gnos_id.split(',')[i]
                    tumor_wgs_alignment_bam_file_name_i = tumor_wgs_alignment_bam_file_name.split(',')[i]
                    #print(tumor_wgs_alignment_bam_file_name_i) 
                    g.write(f'\n{final_index}\t{donor_unique_id}\t{dcc_project_code}\t{submitter_donor_id}\t{icgc_donor_id}\t{tumor_wgs_submitter_specimen_id_i}\t{tumor_wgs_icgc_specimen_id_i}\t{tumor_wgs_submitter_sample_id_i}\t{tumor_wgs_icgc_sample_id_i}\t{tumor_wgs_aliquot_id_i}\t{tumor_wgs_alignment_gnos_repo_i}\t{tumor_wgs_alignment_gnos_id_i}\t{tumor_wgs_alignment_bam_file_name_i}\ttumor\t{tumor_wgs_specimen_count}')
            linecount += 1

