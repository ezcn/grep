import pandas as pd
import argparse

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", help="path fathmm non coding annotated file ",type=str, required= True)
	parser.add_argument("-csq", help="path to csq file ",type=str, required= True)
	parser.add_argument("-o", help="path to output file  ",type=str, required= True)
	args = parser.parse_args()

	noncod={}
	with open(args.f , 'r') as myd:
		firstLine= myd.readline()
		for line in myd:
			elements=line.strip().split("\t")
			mychr=elements[0]; mypos=elements[1]; myalt=elements[3]; myscore=elements[4]
			mykey=mychr + ":" + mypos + ":/" + myalt
			noncod[mykey]=myscore
	df=pd.DataFrame.from_dict(noncod.items())
	df.columns=['index_x', 'fathmmScore']
	df2= pd.read_csv(args.csq,sep="\t")
	dfall = pd.merge(df,df2,on="index_x", how="right").fillna("0")
	dfall.to_csv(args.o, sep="\t", index=False)

if __name__ == "__main__":
	main()
