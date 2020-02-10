#rename file based on first and last row
#python code
import multiprocessing
import os
import sys
import time
import threading
import glob

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

def wc_mp_pool(name, blocksize=65536):
    ranges = get_file_offset_ranges(name, blocksize)

    pool = multiprocessing.Pool(processes=len(ranges))
    pool_outputs = pool.map(get_chunk_line_count, ranges)
    pool.close()
    pool.join()

    return sum(pool_outputs)


print(wc_mp_pool(sys.argv[1]))

#create a list of files
print(">>> create a list of files")
filenames = (glob.glob("SNV_key_ch*.tsv"))
#count the lines in each files
print(">>> count the lines in each files")
len_of_dfs = [len(open(filename).readlines()) for filename in filenames]
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
		new = dataframe.filename_new[0]
		dataframe.loc[:,'filename_new'] = "SNV_key_chr{chr}_{start}-{end}.tsv".format(chr=dataframe.chr[0],start=dataframe.position[0],end=dataframe.position[1])
		print(">> changing name file: {old} --->> {new}".format(old=old,new=new))
		os.rename(old,new)
	else:
		print("ERROR: dataframe not contain unique chr in...QUIT!")
		print(dataframe)
		break
print(">>> changing name: END")

print("QUIT")


def chunks(l, n):
    return [l[i:i+int(n)] for i in range(0, len(l), int(n))]

def job_1(job_id,data_slice,return_dict):
    #are they share somethings?
    return_dict[job_id] = len(open(filename).readlines())

def dispatch_jobs_1(data, job_number=40):
    total = len(data)
    chunk_size = total / job_number
    slice = chunks(data, chunk_size)
    jobs = []
    manager = multiprocessing.Manager()
    return_dict = manager.dict()

    for job_id, data in enumerate(slice):
        j = multiprocessing.Process(target=job_number_1, args=(job_id, data, return_dict))
        jobs.append(j)
    for j in jobs:
        j.start()
    for p in jobs:
        p.join()

    return return_dict
