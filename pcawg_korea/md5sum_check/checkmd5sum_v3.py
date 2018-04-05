import os
import pathlib
import hashlib
import concurrent.futures
import argparse
import re


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--folder', required=True, help='Folder path to check md5sum')
    parser.add_argument('-p', '--suffix_pattern', required=False, help='suffix pattern to search for in the fastq folder', default='*.fastq.gz')

    args = vars(parser.parse_args())

    folder = args['folder']
    suffix_pattern = args['suffix_pattern']

    return folder, suffix_pattern



def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def main():
    folder, suffix_pattern = argument_parser()

    local_md5 = open( 'localmd5.txt', 'w')
    filelist = []

    for root, directories, filenames in os.walk(folder):
        for filename in filenames:
            abspath = os.path.join(root, filename)
            if re.search(string=abspath, pattern=suffix_pattern) != None:
                if not os.path.islink(abspath):
                    if not os.path.isdir(abspath) and os.path.isfile(abspath):
                        #print(abspath)
                        filelist.append(abspath)


    with concurrent.futures.ProcessPoolExecutor() as executor:
        for file, md5value in zip(filelist, executor.map(md5, filelist)):
            local_md5.write(str(md5value) + '\t' + os.path.basename(file) + '\n')

    local_md5.close()

    return 0

if __name__=='__main__':
    main()
