import pandas as pd
import os
def attach_center_write(path, centerName):
    data = pd.read_table(path, names=['md5sum', 'bam'])
    data['centerName'] = centerName
    output = os.path.basename(path) + '.center'
    data.to_csv(output, sep='\t', header=False, index=False)
    return data

attach_center_write('/home/users/cjyoon/Projects/pcawg/received_md5sum/ewha/localmd5_CESC.filtered.txt', 'EWHA')
attach_center_write('/home/users/cjyoon/Projects/pcawg/received_md5sum/ewha/localmd5_UCEC.filtered.txt', 'EWHA')
attach_center_write('/home/users/cjyoon/Projects/pcawg/received_md5sum/snuh/localmd5_snu.txt', 'SNUH')
attach_center_write('/home/users/cjyoon/Projects/pcawg/received_md5sum/snuph/localmd5-BRCA_UK.txt', 'SNUPH')
attach_center_write('/home/users/cjyoon/Projects/pcawg/localmd5.txt', 'KAIST')

