import sqlite3
import pandas as pd
import re, sys, argparse

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-db", help="path to csq db  ",type=str,required=True)
	parser.add_argument("-cg", help="path to control genes file ",type=str,required=True)
	parser.add_argument("-sl", help="list of sampes id ",type=str,required=True)
	parser.add_argument("-pathTodir", help="path to file with sample allele counts directory  ",type=str, required= True)
	parser.add_argument("-chrom", help="chromosome name  ",type=str, required= True)
	parser.add_argument("-ac", help=" allele count >= of   ",type=int , required= True)
	parser.add_argument("-o", help="path to output file  ",type=str, required= True)
	args = parser.parse_args()

	#### access to db
	conn = sqlite3.connect(args.db)
	c = conn.cursor()
	### open file with candidate genes
	df=pd.read_table(args.cg, sep="\t")
	df.columns = df.columns.str.strip()
	df.to_sql("controlGenes", conn)
	###### find candidate genes in csq files
	query = "SELECT * FROM noncodjoin JOIN controlGenes ON noncodjoin.Gene = controlGenes.Gene;"
	dfCtr = pd.read_sql_query(query, conn)

	##### add sample and AltCounts
	listSamples = [line.rstrip('\n') for line in open(args.sl)]
	allcontrol=pd.DataFrame()
	for ss in listSamples:
		tmpPda=dfCtr
		tmpSamp=pd.read_table('%s/%s.%s_counts.tsv' %( args.pathTodir, ss, args.chrom) )
		if "chr" in tmpSamp.loc[0].key:
			tmpSamp.key = tmpSamp.key.str.lstrip("chr")
		if "chr" in tmpPda.loc[0].index_x:
			tmpPda.index_x = tmpPda.index_x.str.lstrip("chr")
		tmpSamp = tmpSamp[ (tmpSamp['ALTcount']>= args.ac )]
		df=tmpPda.reset_index(drop=True).merge(tmpSamp, left_on='index_x', right_on='key').drop("key",axis=1)
		allcontrol=pd.concat([allcontrol, df])	
	allcontrol.to_csv(args.o, sep = "\t", na_rep= "NA", index = False)

if __name__ == "__main__":
	main()