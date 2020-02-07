#rename file based on first and last row
#python code

import glob,os,sys
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
