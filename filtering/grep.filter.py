#!/usr/bin/python3
import pandas as pd
import glob, argparse


def mergeThis(mydir, ID_sample):
	df = pd.DataFrame()
	fileList=glob.glob(mydir)
	for f in fileList:
		tmp = pd.read_csv(f,sep="\t", index_col='index_x')
		tmp = tmp[["csqAllele", "impact","element_id", "type", "csqCount", "gene_id", "gene_symbol","most_severe_consequence", "start", "id", "rare", "commonCSQall", "EmbryoDev", "DDD", "Lethal", "Essential", "Misc", "soScore", "pLI", "chr", "CADD"]]
		df = df.append(tmp)
	df.loc[:,"sample"] = ID_sample
	df.loc[:,"sumGene"] = df["EmbryoDev"]+ df["DDD"]+ df["Lethal"]+ df["Essential"]+ df["Misc"]
	df = df[["sample","csqAllele", "impact","element_id", "type", "csqCount", "gene_id", "gene_symbol","most_severe_consequence", "start", "id", "rare", "commonCSQall", "EmbryoDev", "DDD", "Lethal", "Essential", "Misc", "soScore", "pLI", "chr", "CADD","sumGene"]]
	#df.to_csv("cicci", sep="\t", index=True)
	return df

def main():
        parser = argparse.ArgumentParser()
        parser.add_argument("-f", help="path to  input dir  ",type=str,required=True)
        #parser.add_argument("-so", help="threshold for SOTerm Impact  ", type=float,required=True)
        parser.add_argument("-r", help="threshold for rare variant definition ", type=str,required=True)
        parser.add_argument("-pli", help="threshold for  pLI score  ",type=float, required= True)
        parser.add_argument("-g", help=" number of gene lists  ",type=float , required= True)
        parser.add_argument("-cadd", help=" treshold for CADD score  ",type=float , required= True)
        parser.add_argument("-count", help=" treshold for csq allele count  ",type=float , required= True)
        parser.add_argument("-sample", help="id of the sample  ",type=str, required= True)
        parser.add_argument("-o", help="path to output file  ",type=str, required= True)
        #parser.add_argument("-e", help="path to error file",type=str,required=True)
        args = parser.parse_args()
        #sys.stdout=open(args.o, 'w')   

        #df_allSamples=pd.DataFrame()
        #for ID_sample in samples_grep:
        df= mergeThis(args.f, args.sample)
	
        ##~~~  filter in genic regions
        #df_genic=df.loc[df['type'] =='genic'] 
        df_filtered_g = df[ (df['type']== 'genic') & (df.impact!="MODIFIER") & (df.rare != args.r) & (df.pLI >= args.pli) & ( df.sumGene >= args.g) & (df.CADD >= df.CADD.quantile(args.cadd)) & (df.csqCount >= args.count)]
        #df_allSamples = df_filtered_g.append(df_filtered_g,sort=True)
        df_filtered_g.to_csv(args.o, sep="\t",index=True)

                ##~~~  filter in regulatory regions
                #df_regulatory=df.loc[df['type'] =='regulatory'] 

                ##~~~   filter in intergenic 



if __name__ == "__main__":
        main()
