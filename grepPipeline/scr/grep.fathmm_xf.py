import pandas as pd
import argparse

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-n", help="path to noncoding score file  ",type=str, required= True)
	parser.add_argument("-i", help="path to input file  ",type=str, required= True)
	parser.add_argument("-o", help="path to output file  ",type=str, required= True)
	args = parser.parse_args()

	df1=pd.read_csv(args.n, sep="\t", header=None)
	df2=pd.read_csv(args.i, sep="\t", header=None)
	Allscore=pd.merge(df1,df2)
	Allscore.columns=['chr', 'pos', 'ref', 'alt', 'score']
	Allscore.to_csv(args.o, sep="\t", index=False)

if __name__ == "__main__":
	main()
