#!/usr/bin/python3
import pandas as pd
import glob, argparse


def selectFiles(mydir, listID):
    """ choose files if ID in list id is in the filename""" 
    fileList=[]
    for ind in listID: 
        fileList += [f for f in glob.glob(mydir) if re.search(ind, f) ]
    return fileList

def mergeThat(fileList):
    df = pd.DataFrame()
    for f in fileList:
        tmp = pd.read_table(f) 
        df = df.append(tmp)
    return df

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
        #HGDP arguments
        parser.add_argument("-rf", help="path to  input dir of reference files  ",type=str,required=True)
        parser.add_argument("-l", help="list of id ",type=str,required=True)
        parser.add_argument("-n", help="number of individual to sample ",type=int,required=True)
        parser.add_argument("-i", help="number of iterations ",type=int,required=True) 
        parser.add_argument("-gt", help="threshold for excluding genes ",type=float,required=True) 

        #GREP arguments
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

        #~~~  HGDP 
        listID = [line.rstrip('\n') for line in open(args.l)]
    
        cycle=0
        while cycle < args.i: 
                cycle+=1

                #choose a random sample     
                sampleToConsider=random.sample(listID, args.n)
                #print (sampleToConsider)
    
                # build myd data frame with data from  sample choosen
                fileToConsider=selectFiles(args.rf, sampleToConsider)
                myd = mergeThat(fileToConsider)
                #myd.to_csv("cicci", sep="\t", index=True)
 
                #count each gene only once per sample -> gene symbol, in how many samples it is found 
                tmpg=myd[[ 'gene_symbol', 'sample']].drop_duplicates().groupby(['gene_symbol']).count().transform(lambda x: x /float(args.n) )
                if cycle==1:  genesPerSample=tmpg 
                else: genesPerSample=genesPerSample.join(tmpg, on='gene_symbol', how='outer', lsuffix='_genesPerSample', rsuffix='_tmpg').fillna(0)

        #variantsPerGene.to_csv('ciccivar', sep='\t', index=False)    
        genesPerSample['GrandMean']=genesPerSample.sum(axis=1 )/float(args.i)
        genesToDiscard=genesPerSample[genesPerSample['GrandMean']> float(args.gt)]

        #genesPerSample.to_csv('ciccigene', sep='\t', index=False)
        #genesToDiscard['gene_symbol'].to_csv('ciccigenediscard', sep='\t', index=False)


        #~~~  GREP 
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
