import argparse
import os

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sourceMD5', required=True, help='MD5 values from original source')
    parser.add_argument('-l', '--localMD5', required=True, help='MD5 values from local copy')

    args = vars(parser.parse_args())

    sourceMD5 = args['sourceMD5']
    localMD5 = args['localMD5']

    return sourceMD5, localMD5

def md5dict(path):
    '''create dictionary of file -> md5 value from a given md5 summary'''
    md5dictionary = dict()
    with open(path, 'r') as f:
        for line in f:
            md5value, filepath = line.strip().split()
            md5dictionary.update({os.path.basename(filepath):md5value})

    return md5dictionary

def missing_local(source_dictionary, local_dictionary):
    ''' compares local md5 and source md5 and reports files that are missing locally'''

    missing_list = []
    for filepath in source_dictionary.keys():
        if not filepath in local_dictionary.keys():
            missing_list.append(filepath)

    print(f'###These files are missing locally: {len(missing_list)}')
    printlist(missing_list)

    return missing_list

def printlist(list):
    for element in list:
        print(element)
    return 0

def match_md5(source_dictionary, local_dictionary):
    '''compares local md5 with those in the source'''
    nosource = []
    match = []
    mismatch = []
    for localfile in local_dictionary.keys():
        if localfile in source_dictionary.keys():
            if source_dictionary[localfile] == local_dictionary[localfile]:
                match.append(localfile)
            else:
                mismatch.append(localfile)
        else:
            nosource.append(localfile)

    print(f'\n###no source md5: {len(nosource)} ')
    printlist(nosource)

    print(f'\n###matching files: {len(match)}')
    printlist(match)

    print(f'\n###mismatching files: {len(mismatch)}')
    printlist(mismatch)

    return match, mismatch

def main():
    sourceMD5, localMD5 = argument_parser()
    source_dictionary = md5dict(sourceMD5)
    local_dictionary = md5dict(localMD5)

    missing_list = missing_local(source_dictionary, local_dictionary)

    match_md5(source_dictionary, local_dictionary)

    return 0


if __name__=='__main__':
    main()

