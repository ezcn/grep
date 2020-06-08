import sqlite3
import pandas as pd
import re, sys, argparse

def main():
	parser = argparse.ArgumentParser()
	###~~~~ new sql db name
	parser.add_argument("-db", help="path to input vep table file ",type=str,required=True)
	###~~~ input files
	parser.add_argument("-i", help="path to input vep table file ",type=str,required=True)
	###~~~~ databases
	parser.add_argument("-g", help="path to gene file list ",required=True)
	parser.add_argument("-pli", help="path to table of pLI score  ",type=str, required= True)
	parser.add_argument('-fathmmcoding', help="path to fathmm file",type=str, required= True)
	parser.add_argument('-fathmmnc', help="path to fathmm coding file",type=str, required= True)
	####~~~ threshold for rare variants
	parser.add_argument("-r", help="threshold for rare variant definition ", type=float,required=True)
	###~~~ output file
	parser.add_argument("-o", help="path to output file  ",type=str, required= True)
	args = parser.parse_args()


	conn = sqlite3.connect(args.db)  ### create new sql db
	c = conn.cursor()
	###### open vep table with pandas and create table inside grep.db 
	df=pd.read_table(args.i, sep="\t", index_col= "Uploaded_variation")  ####apro file con pandas na_values="-"
	df.columns = df.columns.str.strip()
	df.to_sql("myTable", conn)
	###### open pLI table and create new table inside grep.db
	df=pd.read_table(args.pli, sep="\t")
	df.columns = df.columns.str.strip()
	df.to_sql("pLItable", conn)
	###### open geneList table and create new table inside grep.db
	df=pd.read_table(args.g, sep="\t")
	df.columns = df.columns.str.strip()
	df.to_sql("geneList", conn)
	##### open fathmm coding table and create new table inside grep.db
	df=pd.read_table(args.fathmmcoding, sep="\t")
	df.columns = df.columns.str.strip()
	df.to_sql("fatmTab", conn)
	###### open fathmm non coding
	df=pd.read_table(args.fathmmnc, sep="\t")
	df.columns = df.columns.str.strip()
	df.to_sql("fatmNCtab", conn)

	
	###### join myTable and pLItab on common column creating new table

	c.execute('CREATE TABLE firstjoin (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,Position_in_cDNA integer,Position_in_CDS integer,Position_in_protein integer,Amino_acid_change text,Codon_change text,Co_located_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CLIN_SIG_ text, pLIscore real);')

	c.execute('INSERT INTO firstjoin SELECT myTable.*, pLItable.pLI FROM myTable LEFT JOIN pLItable ON myTable.Feature = pLItable.transcript;')
	conn.commit()


	##### join genelist

	c.execute('CREATE TABLE genesjoin (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,Position_in_cDNA integer,Position_in_CDS integer,Position_in_protein integer,Amino_acid_change text,Codon_change text,Co_located_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CLIN_SIG_ text, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real);')

	c.execute('INSERT INTO genesjoin SELECT firstjoin.*, geneList.EmbryoDev, geneList.DDD, geneList.Lethal, geneList.Essential, geneList.Misc FROM firstjoin LEFT JOIN geneList ON firstjoin.Gene = geneList.ensID;')
	conn.commit()

	####### join fathmm score
	c.execute('CREATE TABLE codjoin (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,Position_in_cDNA integer,Position_in_CDS integer,Position_in_protein integer,Amino_acid_change text,Codon_change text,Co_located_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CLIN_SIG_ text, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, FathmmCod real);')

	c.execute('INSERT INTO codjoin SELECT genesjoin.*, fatmTab.FathmmCoding FROM genesjoin LEFT JOIN fatmTab ON genesjoin.Location ||\':/\'|| genesjoin.Allele = fatmTab.key;')
	conn.commit()
	###### join fathm non coding

	c.execute('CREATE TABLE noncodjoin (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,Position_in_cDNA integer,Position_in_CDS integer,Position_in_protein integer,Amino_acid_change text,Codon_change text,Co_located_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CLIN_SIG_ text, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, FathmmCod real, FathmmNonCod real);')

	c.execute('INSERT INTO noncodjoin SELECT codjoin.*, fatmNCtab.FathmmNonCod FROM codjoin LEFT JOIN fatmNCtab ON codjoin.Location ||\':/\'|| codjoin.Allele = fatmNCtab.key;')
	conn.commit()

	
	###### retrieve frequencies

	thr = args.r #### threshold for rare variants

	c.execute('CREATE TABLE freqTable (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,Position_in_cDNA integer,Position_in_CDS integer,Position_in_protein integer,Amino_acid_change text,Codon_change text,Co_located_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CLIN_SIG_ text, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, FathmmCod real, FathmmNonCod real,rare text);')

	c.execute("INSERT INTO freqTable SELECT *, CASE WHEN AFR_AF <? AND AMR_AF <? AND ASN_AF <? AND EUR_AF <? AND EAS_AF <? AND SAS_AF <? AND AA_AF <? AND EA_AF <? AND gnomAD_AFR_AF <? AND gnomAD_AMR_AF <? AND gnomAD_ASJ_AF <? AND gnomAD_EAS_AF <? AND gnomAD_FIN_AF <? AND gnomAD_NFE_AF <? AND gnomAD_OTH_AF <? AND gnomAD_SAS_AF <? THEN 'true' WHEN AFR_AF ='-' AND AMR_AF ='-' AND ASN_AF ='-' AND EUR_AF ='-' AND EAS_AF ='-' AND SAS_AF ='-' AND AA_AF ='-' AND EA_AF ='-' AND gnomAD_AF ='-' AND gnomAD_AFR_AF ='-' AND gnomAD_AMR_AF ='-' AND gnomAD_ASJ_AF ='-' AND gnomAD_EAS_AF ='-' AND gnomAD_FIN_AF ='-' AND gnomAD_NFE_AF ='-' AND gnomAD_OTH_AF ='-' AND gnomAD_SAS_AF ='-' THEN 'NOB' ELSE 'false' END FROM noncodjoin;" ,(thr,thr,thr,thr,thr,thr,thr,thr,thr,thr, thr,thr,thr,thr,thr,thr,))

	conn.commit()
	
	######## select column and write a csv file
	
	query  = "SELECT Uploaded_variation, Location ,Allele ,Gene ,Feature ,Feature_type ,Consequence,Position_in_cDNA ,Position_in_CDS,Position_in_protein ,Amino_acid_change ,Codon_change,Co_located_variation ,IMPACT,SYMBOL ,STRAND ,SIFT ,PolyPhen ,EXON ,AF,MAX_AF,CLIN_SIG_ , pLIscore , EmbryoDev , DDD, Lethal, Essential, Misc, FathmmCod , FathmmNonCod,rare FROM freqTable WHERE rare = 'NOB' OR rare = 'true';"

	df_final = pd.read_sql_query(query, conn).fillna(0)
	df_final.to_csv(args.o, sep="\t", index=False)

if __name__ == "__main__":
	main()
