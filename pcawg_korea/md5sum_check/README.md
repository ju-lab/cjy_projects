# md5check
## create a local md5 value table in a given folder

`checkmd5sum.py` will create a two-column table with md5sum value of a file in a given folder with a suffix pattern. For example, if you wanted to get a md5sum value of all fastq.gz files in `/path/to/fastq` folder, simply type

```
python checkmd5sum.py -f /path/to/fastq -p '*.fastq.gz' 
```

will write results into `localmd5.txt`

## compare with source md5 values
`md5compare.py` will compare the md5sum value from source (typically provided) and local copy of the same data. 

```
python md5compare.py -s /source/md5value.txt -l /local/md5value.txt
```

This will report result in four categories
- md5 value is missing from source data
- md5 value is missing from local data 
- md5 value of source and local copy matching
- md5 value of source and local copy mismatching



