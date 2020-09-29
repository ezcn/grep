#!/usr/bin/python3
import pandas as pd

samples=["181755","182151","190116","190742","190810","190923","191080","192128","192369","192682","190798","190423","191460","192127","192824","192361"]
samples_grep=["AS006","AS054","AS064","AS074","AS090","AS094"]
chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
chromosomes_auto = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]
csq_impact = pd.read_csv("/home/adriano/github/grep/filtering/csqimpact.tsv",sep="\t")
dSOTermFineRank=VepRankingInfo("/home/adriano/github/grep/filtering/csqimpact.tsv")
soScore = pd.Series(dSOTermFineRank,name="soScore").to_frame().reset_index()
impact = soScore.merge(csq_impact,left_on="index",right_on="SO term")[["soScore","IMPACT"]] 
###########################################à
#mem
df_big = pd.DataFrame()
df_tmp = pd.DataFrame()
for ID_sample in samples:
	for chrs in chromosomes:
		ID_unique = ID_sample+"/"+ID_sample+"."+chrs+".tsv"
		df = pd.read_csv(ID_unique,sep="\t",index_col="index_x")
		df.loc[:,"sample"] = ID_sample
		df.loc[:,"chr"] = chrs
		df = (df.reset_index().merge(impact)).set_index("index_x")
		df_tmp = df_tmp.append(df,sort=True)
	df_tmp = df_tmp[['sample','chr','id','gene_id','gene_symbol','starting_position','most_severe_consequence','csqCount','csqAllele','score_gene_list','rare','soScore','IMPACT']]
	df_tmp.to_csv(ID_sample+"/"+ID_sample+".whole_genome.tsv",sep="\t")
	df_big = df_big.append(df_tmp)
df_big.to_csv("all_samples.whole_genome.tsv",sep="\t")

########################################
#grep 

df_big = pd.DataFrame()
df_tmp = pd.DataFrame()
for ID_sample in samples_grep:
	for chrs in chromosomes_auto:
		try:
			ID_unique = ID_sample+"/"+ID_sample+"."+chrs+".tsv"
			df = pd.read_csv(ID_unique,sep="\t",index_col="index_x")
			df.loc[:,"sample"] = ID_sample
			df.loc[:,"chr"] = chrs
			df = (df.reset_index().merge(impact)).set_index("index_x")
			df_tmp = df_tmp.append(df,sort=True)
		except Exception as e:
			continue
	df_tmp = df_tmp[['sample','chr','id','gene_id','gene_symbol','starting_position','most_severe_consequence','csqCount','csqAllele','score_gene_list','rare','soScore','IMPACT']]
	df_tmp.to_csv(ID_sample+"/"+ID_sample+".whole_genome.tsv",sep="\t")
	df_big = df_big.append(df_tmp)
df_big.to_csv("all_samples.whole_genome.tsv",sep="\t")

#ANALYSIS
#####################
19 2996512 2996512 
19 2995094 2995094 
19 2976849 2976849 
19 2976547 2976547 
19 2977290 2977290 
19 2987502 2987502 
19 2984349 2984349  

#numer of variant per sample
df_big.groupby("sample")["index_x"].nunique()
181755    3892576
182151    4587795
190116    6064576
190423    1297530
190742    7486583
190798    3496410
190810    6813828
190923    5729251
191080    5032020
191460    2437772
192127    1427580
192128    5613075
192361     438909
192369    3324256
192682    3677940
192824    1250550
df_big.groupby("sample").size().mean()
3910665.6

#% of category per samples
c = df_big.groupby(['sample', 'IMPACT']).size().rename("count")
per_sample_cate = c / df_big.groupby("sample").size()
per_sample_cate.name = "percentage"
per_sample_cate.to_csv("percentage_sample_impact.tsv",sep="\t",header=True) 
#file creato

#filtro
df_big_filtered = df_big[(df_big.score_gene_list > df.score_gene_list.quantile(.99)) & (df_big.soScore > df.soScore.quantile(.99)) & (df_big.rare != False)] 
#filtrato
(TLE6_pre.score_gene_list > TLE6_pre.score_gene_list.quantile(.99)) & (TLE6_pre.soScore > TLE6_pre.soScore.quantile(.99)) & (TLE6_pre.rare != False)
#numero di varianti per campione dopo filtro
df_big_filtered.groupby("sample").size()
sample
181755    36640
182151    37140
190116    37338
190423    12850
190742    37986
190798    19380
190810    35508
190923    31999
191080    29720
191460    11628
192127     8367
192128    27324
192361     2928
192369    20976
192682    20349
192824     6114
#According to Ensembl classification 0.25% of all variants are assumed to have high disruptive impact in the gene product. 
df_big.groupby("IMPACT").size() / len(df_big) * 100       

#Preimplantation development arrest
pre_dev_arrest = df_big[df_big["sample"].isin(pre)] 


pre_dev_arrest_filtered = pre_dev_arrest[(pre_dev_arrest.score_gene_list > pre_dev_arrest.score_gene_list.quantile(.95)) & (pre_dev_arrest.soScore > pre_dev_arrest.soScore.quantile(.99)) & (pre_dev_arrest.rare != False)]
pre_dev_arrest_filtered.to_csv("Preimplantation_deve_arrest_rare_variant_filtered.tsv",sep="\t")






array(['RBBP4', 'EOMES', 'HUS1', 'PTPA', 'MTCH2',
       'TSFM', 'AC025165.3', 'EIF5A', 'PPP5C', 'MAD2L2',
       'RPN1', 'POT1', 'NUP214', 'TP53', 'ABCB7',
       'RBM48', 'MED24', 'SPAG5', 'RPL13A', 'CHEK1',
       'ZFP91', 'ZFP91-CNTF', 'BLM', 'PNKP', 'RBMX',
       'CASP9', 'GON4L', 'REV3L', 'LEO1', 'MRPL28',
       'BIRC5', 'STAG2', 'KPNA4', 'FANCC', 'CIAPIN1',
       'HSCB', 'ORC1', 'MRPS28', 'CDC42BPG', 'PTPMT1',
       'SNRNP35', 'PPM1D', 'LUC7L3', 'SAMD4B', 'COQ2',
       'ASCC3', 'SSSCA1', 'NSUN4', 'NUP54', 'BUD31',
       'NOP53', 'NMD3', 'GET4', 'FAM98B'], dtype=object)


pre_dev_arrest_filtered = pre_dev_arrest[(pre_dev_arrest.score_gene_list > pre_dev_arrest.score_gene_list.quantile(.95)) & (pre_dev_arrest.soScore > pre_dev_arrest.soScore.quantile(.99)) & (pre_dev_arrest.rare != False)]

df_big_nodup_filtered = df_big_nodup[(df_big_nodup.score_gene_list > df_big_nodup.score_gene_list.quantile(.95)) & (df_big_nodup.soScore > df_big_nodup.soScore.quantile(.99)) & (df_big_nodup.rare != False)]
df_big_nodup_filtered.to_csv("grep_filtered.tsv",sep="\t",index=False) 

df_big_nodup_filtered

df_big_nodup_filtered.groupby("sample")["index_x"].nunique()    

uniche_shared = df_big_nodup_filtered.groupby("index_x")["sample"].nunique()
uniche_shared[uniche_shared > 1]


190423    2
191080    1
191460    4
192369

#seleziono prima il tipo di campioni e filtro
pre = ["181755",  "182151",  "190116" , "190742",  "190810",  "190923"  ,"191080" ,"192128","192369","192682"] 
df_fai = df[~df["sample"].isin(pre)]
df_fai_nodup = df_fai.reset_index().drop_duplicates()


#vedere per homozigosita
df_fai_nodup_filtered = df_fai_nodup[(df_fai_nodup.csqCount == 2) & (df_fai_nodup.soScore > df_fai_nodup.soScore.quantile(.99)) & (df_fai_nodup.rare != False)]

shared = df_big_fai_nodup_filtered.groupby("index_x")["sample"].nunique()
shared [shared > 1] 
df_big_fai_nodup_filtered.groupby("sample")["index_x"].nunique()

o 3 geni precedentemente associati con embryo arrest? PLK4,  PADI6 (MIM: 610363; GenBank: NM_207421.4) . e  p66Shc

TLE6 (MIM: 612399) 

gene_search = ["'TLE6'","'PLK4'","'PADI6'","'SHC1'"]
df_arrest = df[df["sample"].isin(pre)]  
for gene in gene_search:
	df_tmp = df_arrest[df_arrest.gene_symbol.str.contains(gene)]
	df_tmp = df_tmp[(df_tmp.rare != False) & (df_tmp.most_severe_consequence != "intron_variant")]
	df_tmp = df_tmp.drop_duplicates()
	print(df_tmp.groupby("sample")["index_x"].unique())
	print(df_tmp.groupby("index_x")["sample"].unique())
	df_tmp.to_excel("gene_{x}_variants.xlsx".format(x=gene),index=False)

