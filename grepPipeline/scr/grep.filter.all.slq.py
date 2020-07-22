import sqlite3
import pandas as pd
import re, sys, argparse

def main():
	parser = argparse.ArgumentParser()
	###~~~~ new sql db name
	parser.add_argument("-db", help="path to new db ",type=str,required=True)
	###~~~ input files
	parser.add_argument("-scsq", help="path to input grep csq file ",type=str,required=True)
	parser.add_argument("-ccsq", help="path to input hgdp csq file ",type=str,required=True)
	parser.add_argument("-cl", help="list of control samples id ",type=str,required=True)
	parser.add_argument("-sl", help="list of sampes id ",type=str,required=True)
	###~~~ filtering arguments
	parser.add_argument("-f", help="threshold for rare frequency definition", type=float,required=True)
	parser.add_argument("-type", help=" feature_type(genic , intergenic, regularoty) ", type=str,required=True)
	parser.add_argument("-r", help=" rare variant not equal to", type=str,required=True)
	parser.add_argument("-pli", help="threshold for  pLI score  ",type=float, required= True)
	parser.add_argument("-cadd", help=" treshold for CADD score  ",type=float , required= True)
	parser.add_argument("-g", help=" number of gene lists  ",type=float , required= True)
	###~~~ HGDP arguments
	parser.add_argument("-n", help="number of individual to sample ",type=int,required=True)
	parser.add_argument("-i", help="number of iterations ",type=int,required=True)
	parser.add_argument("-pathTodirCtrl", help="path to file directory  ",type=str, required= True)
	parser.add_argument("-ctgn", help="path to control genes file",type=str,required=True)
	parser.add_argument("-gtd", help="path to genes to discard file",type=str,required=True)

	###~~~ Grep arguments
	parser.add_argument("-pathTodir", help="path to file directory  ",type=str, required= True)
	parser.add_argument("-chrom", help="chromosome name  ",type=str, required= True)
	parser.add_argument("-ac", help=" allele count >= of   ",type=int , required= True)
	#parser.add_argument("-ctrlgen", help="path to hgdp genes to discard file ",type=str,required=True)
	parser.add_argument("-gt", help="threshold for excluding genes ",type=float,required=True)
	parser.add_argument("-maxv", help=" maximum variants per gene ",type=int , required= True)
	parser.add_argument("-o", help="path to output file  ",type=str, required= True)
	args = parser.parse_args()


	conn = sqlite3.connect(args.db)  ### create new sql db
	c = conn.cursor()
	###### open vep table with pandas and create table inside grep.db 
	df_grep=pd.read_table(args.scsq, sep="\t", index_col= "Uploaded_variation")  ####apro file con pandas na_values="-"
	df_grep.columns = df.columns.str.strip()
	df_grep.to_sql("myTable", conn)
	df_grep=pd.read_table(args.ccsq, sep="\t", index_col= "Uploaded_variation")
	df_grep.columns = df.columns.str.strip()
	df_grep.to_sql("myCtr", conn)

	###### retrieve frequencies for control samples
	thr = args.ff

	c.execute('CREATE TABLE freqTable1Ctr (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, FathmmCod real, FathmmNonCod real, rare01 real);')

	c.execute("INSERT INTO freqTable1Ctr SELECT *, CASE WHEN AFR_AF ='-' AND AMR_AF ='-' AND ASN_AF ='-' AND EUR_AF ='-' AND EAS_AF ='-' AND SAS_AF ='-' AND AA_AF ='-' AND EA_AF ='-' AND gnomAD_AF ='-' AND gnomAD_AFR_AF ='-' AND gnomAD_AMR_AF ='-' AND gnomAD_ASJ_AF ='-' AND gnomAD_EAS_AF ='-' AND gnomAD_FIN_AF ='-' AND gnomAD_NFE_AF ='-' AND gnomAD_OTH_AF ='-' AND gnomAD_SAS_AF ='-' THEN 'NOB' WHEN AFR_AF <=? AND AMR_AF <=? AND ASN_AF <=? AND EUR_AF <=? AND EAS_AF <=? AND SAS_AF <=? AND AA_AF <=? AND EA_AF <=? AND gnomAD_AFR_AF <=? AND gnomAD_AMR_AF <=? AND gnomAD_ASJ_AF <=? AND gnomAD_EAS_AF <=? AND gnomAD_FIN_AF <=? AND gnomAD_NFE_AF <=? AND gnomAD_OTH_AF <=? AND gnomAD_SAS_AF <=? THEN 'true' ELSE 'false' END FROM myCtr;" ,(thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,))

	conn.commit()
	
	thr = args.f

	c.execute('CREATE TABLE freqTableCtr (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, FathmmCod real, FathmmNonCod real, rare01 real, rare05 real);')

	c.execute("INSERT INTO freqTableCtr SELECT *, CASE WHEN AFR_AF ='-' AND AMR_AF ='-' AND ASN_AF ='-' AND EUR_AF ='-' AND EAS_AF ='-' AND SAS_AF ='-' AND AA_AF ='-' AND EA_AF ='-' AND gnomAD_AF ='-' AND gnomAD_AFR_AF ='-' AND gnomAD_AMR_AF ='-' AND gnomAD_ASJ_AF ='-' AND gnomAD_EAS_AF ='-' AND gnomAD_FIN_AF ='-' AND gnomAD_NFE_AF ='-' AND gnomAD_OTH_AF ='-' AND gnomAD_SAS_AF ='-' THEN 'NOB' WHEN AFR_AF <=? AND AMR_AF <=? AND ASN_AF <=? AND EUR_AF <=? AND EAS_AF <=? AND SAS_AF <=? AND AA_AF <=? AND EA_AF <=? AND gnomAD_AFR_AF <=? AND gnomAD_AMR_AF <=? AND gnomAD_ASJ_AF <=? AND gnomAD_EAS_AF <=? AND gnomAD_FIN_AF <=? AND gnomAD_NFE_AF <=? AND gnomAD_OTH_AF <=? AND gnomAD_SAS_AF <=? THEN 'true' ELSE 'false' END FROM freqTable1Ctr;" ,(thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,))

	conn.commit()
	###### create sumGene column in control samples df
	c.execute('CREATE TABLE sumgeneTabCtr (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, FathmmCod real, FathmmNonCod real, rare real, sumGene real);')

	c.execute("INSERT INTO sumgeneTabCtr SELECT Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,SYMBOL,STRAND,SIFT,PolyPhen,EXON,MAX_AF,CADD_RAW, CADD_PHRED, pLIscore, EmbryoDev, DDD, Lethal, Essential, Misc, index_x text, FathmmCod, FathmmNonCod,rare, (EmbryoDev + DDD + Lethal + Essential + Misc) FROM freqTableCtr;")

	conn.commit()

	###### create CADD percentile column in control samples df
	c.execute("CREATE TABLE finalTableCtr (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, FathmmCod real, FathmmNonCod real, rare real, sumGene real, caddPercent real);")


	c.execute("INSERT INTO finalTableCtr SELECT * , ROUND(PERCENT_RANK() OVER (ORDER BY CADD_RAW),3) FROM sumgeneTabCtr WHERE CADD_RAW != '-';")

	conn.commit()

	######## filter control samples 
	typ = args.type ; rareThresh = args.r ; pliscore = args.pli ; caddscore = args.cadd ; numgene = args.g

	#c.execute("CREATE TABLE filtro AS SELECT * FROM finalTable WHERE IMPACT != 'MODIFIER' AND feature_type = ? AND rare != ? AND (pLIscore >= ? AND caddPercent >= ? OR sumGene >= ?);" , (typ,rareThresh, pliscore, caddscore, numgene,))
	query = "SELECT * FROM finalTableCtr WHERE IMPACT != 'MODIFIER' AND feature_type = ? AND rare != ? AND (pLIscore >= ? AND caddPercent >= ? OR sumGene >= ?); "
	df_finalCtr = pd.read_sql_query(query,conn, params = (typ,rareThresh, pliscore, caddscore, numgene))
	#query = "SELECT * FROM filtro;"
	#df_finalHg.to_csv()
	#df_finalHg = pd.read_sql_query(query,conn)

	##~~ repeat args.i times on args.n samples from controls
	listControl = [line.rstrip('\n') for line in open(args.cl)]
	cycle=0
	while cycle < args.i:
		cycle+=1
		##~~ choose a random sample form controls of size args.n
		randomSample=random.sample(listControl, args.n)
		##~~ integrate with samples ID and csqAlele count
		control_filtered_allsamples=pd.DataFrame()
		for ss in randomSample:
			tmpdf=df_finalCtr
			tmpSamp=pd.read_table('%s/%s.%s_counts.tsv' %( args.pathTodirCtrl, ss, args.chrom) )
			tmpSamp = tmpSamp[ (tmpSamp['ALTcount']>=args.ac )]
			tmpSamp["sample"] = ss
			df=tmpdf.reset_index(drop=True).merge(tmpSamp, left_on='index_x', right_on='key').drop("key",axis=1)
			control_filtered_allsamples=pd.concat([control_filtered_allsamples, df])

			tmpg=control_filtered_allsamples[[ 'SYMBOL', 'ID']].drop_duplicates().groupby('SYMBOL').count().transform(lambda x: x /float(args.n) )
			if cycle==1:  genesPerSample=tmpg  # create genesPerSample at first cycle
			else: genesPerSample=genesPerSample.join(tmpg, on='SYMBOL', how='outer', lsuffix='_genesPerSample', rsuffix='_tmpg').fillna(0)
		##~~ make GrandMean over args.i and discard genes that on average shows up in args.gt individuals over args.i iterations	
		genesPerSample['GrandMean']=genesPerSample.sum(axis=1 )/float(args.i)
		genesToDiscardControl=genesPerSample[genesPerSample['GrandMean']> float(args.gt)]

		genesPerSample.to_csv(args.ctgn, sep='\t')
		genesToDiscardControl.to_csv(args.gtd, sep='\t')
	
	###### retrieve frequencies Grep samples
	thr = args.ff

	c.execute('CREATE TABLE freqTable1 (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, FathmmCod real, FathmmNonCod real, rare01 real);')

	c.execute("INSERT INTO freqTable SELECT *, CASE WHEN AFR_AF ='-' AND AMR_AF ='-' AND ASN_AF ='-' AND EUR_AF ='-' AND EAS_AF ='-' AND SAS_AF ='-' AND AA_AF ='-' AND EA_AF ='-' AND gnomAD_AF ='-' AND gnomAD_AFR_AF ='-' AND gnomAD_AMR_AF ='-' AND gnomAD_ASJ_AF ='-' AND gnomAD_EAS_AF ='-' AND gnomAD_FIN_AF ='-' AND gnomAD_NFE_AF ='-' AND gnomAD_OTH_AF ='-' AND gnomAD_SAS_AF ='-' THEN 'NOB' WHEN AFR_AF <=? AND AMR_AF <=? AND ASN_AF <=? AND EUR_AF <=? AND EAS_AF <=? AND SAS_AF <=? AND AA_AF <=? AND EA_AF <=? AND gnomAD_AFR_AF <=? AND gnomAD_AMR_AF <=? AND gnomAD_ASJ_AF <=? AND gnomAD_EAS_AF <=? AND gnomAD_FIN_AF <=? AND gnomAD_NFE_AF <=? AND gnomAD_OTH_AF <=? AND gnomAD_SAS_AF <=? THEN 'true' ELSE 'false' END FROM myTable;" ,(thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,))

	conn.commit()

	thr = args.f

	c.execute('CREATE TABLE freqTable (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, FathmmCod real, FathmmNonCod real, rare01 real, rare05 real);')

	c.execute("INSERT INTO freqTable SELECT *, CASE WHEN AFR_AF ='-' AND AMR_AF ='-' AND ASN_AF ='-' AND EUR_AF ='-' AND EAS_AF ='-' AND SAS_AF ='-' AND AA_AF ='-' AND EA_AF ='-' AND gnomAD_AF ='-' AND gnomAD_AFR_AF ='-' AND gnomAD_AMR_AF ='-' AND gnomAD_ASJ_AF ='-' AND gnomAD_EAS_AF ='-' AND gnomAD_FIN_AF ='-' AND gnomAD_NFE_AF ='-' AND gnomAD_OTH_AF ='-' AND gnomAD_SAS_AF ='-' THEN 'NOB' WHEN AFR_AF <=? AND AMR_AF <=? AND ASN_AF <=? AND EUR_AF <=? AND EAS_AF <=? AND SAS_AF <=? AND AA_AF <=? AND EA_AF <=? AND gnomAD_AFR_AF <=? AND gnomAD_AMR_AF <=? AND gnomAD_ASJ_AF <=? AND gnomAD_EAS_AF <=? AND gnomAD_FIN_AF <=? AND gnomAD_NFE_AF <=? AND gnomAD_OTH_AF <=? AND gnomAD_SAS_AF <=? THEN 'true' ELSE 'false' END FROM freqTable1;" ,(thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,))

	conn.commit()
	###### create sumGene column grep samples
	c.execute('CREATE TABLE sumgeneTab (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, FathmmCod real, FathmmNonCod real, rare real, sumGene real);')

	c.execute("INSERT INTO sumgeneTab SELECT Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,SYMBOL,STRAND,SIFT,PolyPhen,EXON,MAX_AF,CADD_RAW, CADD_PHRED, pLIscore, EmbryoDev, DDD, Lethal, Essential, Misc, index_x text, FathmmCod, FathmmNonCod,rare, (EmbryoDev + DDD + Lethal + Essential + Misc) FROM freqTable;")

	conn.commit()

	###### create CADD percentile column grep samples
	c.execute("CREATE TABLE finalTable (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, FathmmCod real, FathmmNonCod real, rare real, sumGene real, caddPercent real);")


	c.execute("INSERT INTO finalTable SELECT * , ROUND(PERCENT_RANK() OVER (ORDER BY CADD_RAW),3) FROM sumgeneTab WHERE CADD_RAW != '-';")

	conn.commit()

	########## FILTERING grep samples

	typ = args.type ; rareThresh = args.r ; pliscore = args.pli ; caddscore = args.cadd ; numgene = args.g

	#c.execute("CREATE TABLE filtro AS SELECT * FROM finalTable WHERE IMPACT != 'MODIFIER' AND feature_type = ? AND rare != ? AND (pLIscore >= ? AND caddPercent >= ? OR sumGene >= ?);" , (typ,rareThresh, pliscore, caddscore, numgene,))
	query = "SELECT * FROM finalTable WHERE IMPACT != 'MODIFIER' AND feature_type = ? AND rare != ? AND (pLIscore >= ? AND caddPercent >= ? OR sumGene >= ?); "
	df_final = pd.read_sql_query(query,conn, params = (typ,rareThresh, pliscore, caddscore, numgene))
	#query = "SELECT * FROM filtro;"
	df_final.to_csv()
	df_final = pd.read_sql_query(query,conn)
	
	listSamples = [line.rstrip('\n') for line in open(args.sl)]
	ssall=pd.DataFrame()
	for ss in listSamples:
		tmpgrep=df_final
		tmpSamp=pd.read_table('%s/%s.%s_counts.tsv' %( args.pathTodir, ss, args.chrom) )
		tmpSamp = tmpSamp[ (tmpSamp['ALTcount']>= args.ac )]
		df=tmpgrep.reset_index(drop=True).merge(tmpSamp, left_on='index_x', right_on='key').drop("key",axis=1)
		ssall=pd.concat([ssall, df])
		#df_final.to_csv("filtropdaslq" , sep="\t", index=False)
		

	ssall.to_sql("grepFilter", conn)
	#df=pd.read_table(args.ctrlgen, sep="\t")
	#df.columns = df.columns.str.strip()
	genesPerSample.to_sql("genesMean", conn)

	hgdpMean = args.gt
	c.execute("CREATE TABLE noCtrlGenes AS SELECT grepFilter.* FROM grepFilter LEFT JOIN genesMean ON grepFilter.SYMBOL = genesMean.gene_symbol WHERE GrandMean < ?; " , (hgdpMean,))

	variants = args.maxv
	c.execute("CREATE TABLE variantsPerGene AS SELECT index_x, SYMBOL ,COUNT (DISTINCT index_x) FROM noCtrlGenes GROUP BY SYMBOL HAVING COUNT (DISTINCT index_x) >=?;" ,(variants,))

	c.execute("DELETE FROM noCtrlGenes WHERE EXISTS (SELECT * FROM variantsPerGene WHERE variantsPerGene.SYMBOL = noCtrlGenes.SYMBOL);")
	conn.commit()
	query = "SELECT * FROM noCtrlGenes;"
	df_end = pd.read_sql_query(query,conn)
	df_end.to_csv(args.o, sep = "\t", index = False)


if __name__ == "__main__":
	main()