""" since naive searching returns some duplicated symbolic links, removed those for analsysis"""
import os 
pcawg_normal_paths = '/nfs/users/nfs_y/ysj/lustreScratch112/cjy_mt_vaf/pcawg_normal_paths.txt'
pcawg_normal_paths_list = []

pcawg_tumour_paths = '/nfs/users/nfs_y/ysj/lustreScratch112/cjy_mt_vaf/pcawg_tumour_paths.txt'
pcawg_tumour_paths_list = []


with open(pcawg_normal_paths, 'r') as f:
    for line in f:
        normalpaths = line.strip()
        pcawg_normal_paths_list.append(os.path.realpath(normalpaths))


with open(pcawg_tumour_paths, 'r') as f:
    for line in f:
        tumourpaths = line.strip()
        pcawg_tumour_paths_list.append(os.path.realpath(tumourpaths))


with open('pcawg_normal_abspath_removed_duplicates.txt', 'w') as f:
    for path in set(pcawg_normal_paths_list):
        f.write(path + '\n')

with open('pcawg_tumour_abspath_removed_duplicates.txt', 'w') as f:
    for path in set(pcawg_tumour_paths_list):
        f.write(path + '\n')


