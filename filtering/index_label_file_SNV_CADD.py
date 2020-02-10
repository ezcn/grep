#rename file based on first and last row
#python code
import multiprocessing
import os
import sys
import time
import threading
import glob
import concurrent.futures as futures
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import numpy as np


def get_chunk_line_count(ranges):
    name, start, stop, blocksize = ranges
    left = stop - start

    def blocks(f, left):
        while left > 0:
            b = f.read(min(left, blocksize))
            if b:
                yield b
            else:
                break
            left -= len(b)

    with open(name, 'r') as f:
        f.seek(start)
        return sum(bl.count('\n') for bl in blocks(f, left))

def get_file_offset_ranges(name, blocksize=65536, m=1):
    fsize = os.stat(name).st_size
    chunksize = (fsize // multiprocessing.cpu_count()) * m
    n = fsize // chunksize

    ranges = []
    for i in range(0, n * chunksize, chunksize):
        ranges.append((name, i, i + chunksize, blocksize))
    if fsize % chunksize != 0:
        ranges.append((name, ranges[-1][2], fsize, blocksize))

    return ranges

def wc_proc_pool_exec(name, blocksize=65536):
    ranges = get_file_offset_ranges(name, blocksize)

    with ProcessPoolExecutor(max_workers=len(ranges)) as executor:
        results = [executor.submit(get_chunk_line_count, param) for param in ranges]
        return sum([future.result() for future in futures.as_completed(results)])

#"equal to wc -l"
#print(wc_proc_pool_exec(sys.argv[1]))

#create a list of files
print(">>> create a list of files")
filenames = (glob.glob("SNV_key_ch*.tsv"))
#count the lines in each files
print(">>> count the lines in each files")
###### sequencial format
#len_of_dfs = [len(open(filename).readlines()) for filename in filenames]
###### parallel format
len_of_dfs = [wc_proc_pool_exec(filename) for filename in filenames]
#read only first and last row in each file
print(">>> read only first and last row in each file")
list_of_dfs = [pd.read_csv(filename,skiprows=range(2,len_of_dfs[n]-1), header=0,sep="\t") for n,filename in enumerate(filenames)]
#change name
print(">>> changing name: START")
for dataframe, filename in zip(list_of_dfs, filenames):
	dataframe.loc[:,'filename_old'] = filename
	dataframe.loc[:,"position"] = (dataframe.key.str.split(":",expand=True))[1]
	dataframe.loc[:,"chr"] = (dataframe.key.str.split(":",expand=True))[0]
	if dataframe.chr.nunique() == 1:
		old = dataframe.filename_old[0]
		dataframe.loc[:,'filename_new'] = "range_chr{chr}_{start}-{end}.tsv".format(chr=dataframe.chr[0],start=dataframe.position[0],end=dataframe.position[1])
		new = dataframe.filename_new[0]
		print(">> changing name file: {old} --->> {new}".format(old=old,new=new))
		os.rename(old,new)
	else:
		print("ERROR: dataframe not contain unique chr in...QUIT!")
		print(dataframe)
		break
print(">>> changing name: END")


print(">>> Creating indix...")
#create dataframe info / indice
filenames_new = pd.Series(glob.glob("range_chr*.tsv"),name="file_name").to_frame()
ranges = filenames_new.file_name.str.split("_",expand=True)[2].str.split(".",expand=True)[0].str.split("-",expand=True)
filenames_new.loc[:,"chr"] = filenames_new.file_name.str.split("_",expand=True)[1]
filenames_new.loc[:,"lows"] = (ranges[0]).astype(int)
filenames_new.loc[:,"ups"] = (ranges[1]).astype(int)
lows = filenames_new["lows"].values  # the lower bounds
ups = filenames_new["ups"].values # the upper bounds

print("QUIT")

#### example 
'''
chrc = "chr1"
x = 6773748
filenames_new[(filenames_new["chr"] == chrc) & (filenames_new["lows"] <= x) & (filenames_new["ups"] >= x)]
'''





