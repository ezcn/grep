#!/usr/bin/python3
import pandas as pd
import numpy as np
import os,argparse
from tqdm import tqdm 
from glob import iglob  
import importlib.util,sys
from bisect import bisect_left
from collections.abc import MutableMapping

class RangeBisection(MutableMapping):
    """Map ranges to values

    Lookups are done in O(logN) time. There are no limits set on the upper or
    lower bounds of the ranges, but ranges must not overlap.

    """
    def __init__(self, map=None):
        self._upper = []
        self._lower = []
        self._values = []
        if map is not None:
            self.update(map)

    def __len__(self):
        return len(self._values)

    def __getitem__(self, point_or_range):
        if isinstance(point_or_range, tuple):
            low, high = point_or_range
            i = bisect_left(self._upper, high)
            point = low
        else:
            point = point_or_range
            i = bisect_left(self._upper, point)
        if i >= len(self._values) or self._lower[i] > point:
            raise IndexError(point_or_range)
        return self._values[i]

    def __setitem__(self, r, value):
        lower, upper = r
        i = bisect_left(self._upper, upper)
        if i < len(self._values) and self._lower[i] < upper:
            raise IndexError('No overlaps permitted')
        self._upper.insert(i, upper)
        self._lower.insert(i, lower)
        self._values.insert(i, value)

    def __delitem__(self, r):
        lower, upper = r
        i = bisect_left(self._upper, upper)
        if self._upper[i] != upper or self._lower[i] != lower:
            raise IndexError('Range not in map')
        del self._upper[i]
        del self._lower[i]
        del self._values[i]

    def __iter__(self):
        yield from zip(self._lower, self._upper)

def indexing(index_file):
    index_dict = {}
    for i,k in index_file.iterrows():
        chro=k["chr"]
        file_name=k["file_name"]
        lows=k["lows"]
        ups=k["ups"]
        if chro not in index_dict:
            index_dict[chro] = RangeBisection({(lows,ups) : file_name})
        else:
            index_dict[chro][(lows,ups)] = file_name
    return index_dict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", help="path to  input dir samples",type=str,required=True)
    parser.add_argument("-index", help="path to  index  file ",type=str,required=True)
    parser.add_argument("-cadd", help="path to CADD file list ",required=True)
    parser.add_argument("-chr", help="specify chromosomes, [only numbers] or X or Y",required=True)
    #parser.add_argument("-output", help="path to output file  ",type=str, required= True)
    args = parser.parse_args()

    index_file = args.index     #"/data/resources/CADD_index/index_file_CADD.tsv"
    input_sample = args.input     #"/data/research/NGS/results/grep"
    CADD_files = args.cadd     #"/home/data/resources/CADD_index"
    #out_file = args.output    #"/data/research/NGS/results/grep/"
    chro = args.chr #10

    print("")
    print("Input FILEs : ",input_sample)
    print("index FILE : ",index_file)
    print("CADD FILEs : ",CADD_files)
    print("chr selected : chr",chro)

    index_file = pd.read_csv(index_file,sep="\t") 
    index_dict = indexing(index_file)
    
    if os.path.isfile(input_sample):
        rootdir_glob = (os.path.dirname(input_sample)+"/*chr{chr}.tsv".format(chr=chro)).replace("//","/")
        file_list = [f for f in iglob(rootdir_glob, recursive=True) if os.path.isfile(f)]
    elif os.path.isdir(input_sample):
        rootdir_glob = (input_sample+"/**/*chr{chr}.tsv".format(chr=chro)).replace("//","/")
        file_list = [f for f in iglob(rootdir_glob, recursive=True) if os.path.isfile(f)]

    #print(os.getcwd())
    sys.path.append(CADD_files)
    os.chdir(CADD_files)

    for pos,f in enumerate(file_list):
        print("processing file >>>  "+f)
        df_final = pd.read_csv(f,sep="\t",low_memory=False)
        df_final.drop_duplicates(inplace=True)
        #get info index
        info_key = df_final["index_x"].str.split(":",expand=True)
        df_final.loc[:,"chr"] = info_key[0]
        df_final.loc[:,"position"] = info_key[1]
        #load info index file insert PATH in /lustrehome/enza/ezcngit/
        #index_file = pd.read_csv("/lustrehome/enza/ezcngit/index_file_CADD.tsv",sep="\t")
        print(">>>> get CADD file index ")
        #get the name of cadd file to use
        def get_CADD_file(x):
            try:
                return index_dict["chr"+x["chr"]][int(x["position"])]
            except Exception as e:
                return np.nan

        def get_CADDscore(df):
            try:
                f = df["file_index"].replace(".py","")
                d = __import__(f).d
                return d[df["index_x"]]
            except Exception as e:
                return np.nan

        df_final.loc[:,"file_index"] = df_final.apply(get_CADD_file,axis=1)
        print(">>>> get CADD score ")
        #get CADD score FINALLY MUCH FASTER with __pycache__!
        to_import = df_final.file_index.dropna().unique()
        how_many = list(range(len(to_import)))
        if pos == 0:
            print(">>>> importing index ")
            with tqdm(total=len(how_many),ncols=80) as pbar:
                for cache in to_import:
                    __import__(cache.replace(".py",""))
                    pbar.update(1)
        df_final.loc[:,"CADD"] = df_final.apply(get_CADDscore,axis=1)
        df_final.drop(["position","file_index"],axis=1,inplace=True)
        df_final.to_csv(f,sep="\t",index=False)
        print(">>>> output file: ",f)
        print()

if __name__ == "__main__":
    main()

