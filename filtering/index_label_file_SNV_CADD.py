


for ogni_elemento in directory:
	df = pd.read_csv(ogni_elemento,sep="\t",nrows=1)


count= 
df=pd.read_csv("SNV_key_ch4207.tsv", skiprows=range(2,count-1), header=0,sep="\t")

import glob
filenames = (glob.glob("SNV_key_ch*.tsv"))[0:10]
len_of_dfs = [len(open(filename).readlines()) for filename in filenames]
list_of_dfs = [pd.read_csv(filename,skiprows=range(2,len_of_dfs[n]-1), header=0,sep="\t") for n,filename in enumerate(filenames)]
for dataframe, filename in zip(list_of_dfs, filenames):
	dataframe.loc[:,"position"] = (dataframe.key.str.split(":",expand=True))[1]
	dataframe.loc[:,"chr"] = (dataframe.key.str.split(":",expand=True))[0]
	if dataframe.chr.nunique() == 1:
		dataframe.loc[:,'filename'] = "SNV_key_chr{chr}_{start}-{end}.tsv".format(chr=dataframe.chr[0],start=dataframe.position[0],end=dataframe.position[1])

