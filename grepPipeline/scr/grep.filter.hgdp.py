#!/usr/bin/python3
import pandas as pd
import glob, argparse,  subprocess



def takeInfoFromVcf(vcf, key , sampleToConsider):
    #TO DO: check concordance ref alt and overlapping multiple positions 
    (chrom,pos,alternate) = key.split(":")
    cmd = "tabix -h  %s %s:%d-%d | tail -2" % (vcf, chrom, int(pos), int(pos))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    result, err = proc.communicate()
    if err: raise IOError("** Error running %s key for %s on %s" % (keyString, db))
    vcfout =[x.rstrip().split('\t') for  x in result.decode().split('\n')]
    column2retain=[vcfout[0].index(ind) for ind in sampleToConsider]
    genotypesToConsider=[vcfout[1][indx].split(":")[0] for indx in column2retain] 
    refAllele=vcfout[0][3];  altAlleles=vcfout[0][4]
    return genotypesToConsider, refAllele,  altAlleles
 
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def Freq_CSQ_REF_ALT (csqAllele, refAllele, altAlleles, missing_data_format, genotypeslist):
    """csqAllele = type: string, consequence allele  from vep  
    refAllele = type: string, reference allele from variant calling
    altAlleles = type: string, comma-separated list of alternate allele from variant calling
    nbAploidSamples : calculated 
    GTfields = type: list, list of genotypes ["0/0", "0/1", ".", "./."]
    hetGenotypes = type: int, heterozygosity in samples"""
    myres=False
    listAlt=altAlleles.split(",")   
    listAll=[refAllele]+ listAlt
    stringOfGenotypes=""; nbHaploidSamples=0
    for item in (genotypeslist):    
        if item != missing_data_format: 
            stringOfGenotypes+=item; nbHaploidSamples+=2
    CountAlleles=[]
    for i in range(len(listAll)):  # 0 for REF, 1 for ALT1, 2 for ALT2 ...
        CountAlleles.append(stringOfGenotypes.count(str(i)))
    dAllele=dict(zip(listAll,CountAlleles))
    if nbHaploidSamples!=0:
        freqREF="{0:4.2f}".format(dAllele[refAllele]/nbHaploidSamples)
        freqAlt=[]
        for i in listAlt: 
            freqAltTemp="{0:4.2f}".format(dAllele[i]/nbHaploidSamples)
            freqAlt.append(freqAltTemp)
        #~~~ CSQ 
        if csqAllele in dAllele:
            csqAllCount=dAllele[csqAllele]
            freqCsq="{0:4.2f}".format(csqAllCount/nbHaploidSamples) 
        else: freqCsq='NA'
        myres= [freqCsq, freqREF, freqAlt]
    return myres

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def selectFiles(mydir, listID):
    """ choose files if ID in list id is in the filename""" 
    fileList=[]
    for ind in listID: 
        fileList += [f for f in glob.glob(mydir) if re.search(ind, f) ]
    return fileList

#def mergeThat(fileList):
    

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def filter (df, genometype, thresold_rare, threshold_pli, threshold_sumgene, threshold_cadd): 
    df_filtered = df[ (df['type']== genometype) & (df.impact!="MODIFIER") & (df.rare != thresold_rare) & (df.pLI >= threshold_pli) & ( df.sumGene >= threshold_sumgene) & (df.CADD >= df.CADD.quantile(threshold_cadd))] #& (df.csqCount >= args.count)]
    return df_filtered

#FUNCTIONS END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



def main():
        parser = argparse.ArgumentParser()
        #HGDP arguments
        parser.add_argument("-ccsq", help="path to control csq  ",type=str,required=True)
        parser.add_argument("-rvcf", help="path to reference vcf  ",type=str,required=True)
        parser.add_argument("-cl", help="list of reference id ",type=str,required=True)
        parser.add_argument("-n", help="number of individual to sample ",type=int,required=True)
        parser.add_argument("-i", help="number of iterations ",type=int,required=True) 
        parser.add_argument("-gt", help="threshold for excluding genes ",type=float,required=True) 

        #GREP arguments
        parser.add_argument("-scsq", help="path to control csq  ",type=str,required=True)
        parser.add_argument("-svcf", help="path to reference vcf  ",type=str,required=True)
        parser.add_argument("-sl", help="list of reference id ",type=str,required=True)

        #FILTERING arguments    
        parser.add_argument("-r", help="threshold for rare variant definition ", type=str,required=True)
        parser.add_argument("-pli", help="threshold for  pLI score  ",type=float, required= True)
        parser.add_argument("-g", help=" number of gene lists  ",type=float , required= True)
        parser.add_argument("-cadd", help=" treshold for CADD score  ",type=float , required= True)
        #parser.add_argument("-count", help=" treshold for csq allele count  ",type=float , required= True)

        parser.add_argument("-o", help="path to output file  ",type=str, required= True)
        parser.add_argument("-e", help="path to error file",type=str,required=True)
        args = parser.parse_args()
        #sys.stdout=open(args.o, 'w')   

        #~~~  HGDP 
        #Random sampling args.i times of args.n individual to figure out genes that show up on average args.gt% times over args.i  iterations. Annotate genesToDiscard in a list and exclude these genes from results 
        control = pd.read_table(args.rcsq) 
        control_filtered=filter(control, 'genic', args.r, args.pli, args.g, args.cadd)
        
        listControl = [line.rstrip('\n') for line in open(args.cl)]
        cycle=0
        while cycle < args.i: 
                cycle+=1

                #choose a random sample     
                sampleToConsider=random.sample(listControl, args.n)
                #print (sampleToConsider)

                ####################################################
                #integrate with allele freqency in sample to consider
                mapperAF = {}
                for key in control_filtered.index_x.unique():
                    (chrom,pos,alternate) = key.split(":")
                    genotypesToConsider, refAllele,  altAlleles = makegenotypelist(args.rvcf, key , sampleToConsider)
                    mapperAF.update(Freq_CSQ_REF_ALT (csqAllele, refAllele, altAlleles, missing_data_format, genotypeslist))    
                df_info["AF"] = control_filtered.index_x.map(mapperAF)

                #count each gene only once per sample -> gene symbol, in how many samples it is found 
                #subset by samples in sampleToConsider, subset by loci with AF>0
                tmpg=control_filtered.filter('AF'>0)[[ 'gene_symbol', 'sample']].drop_duplicates().groupby(['gene_symbol']).count().transform(lambda x: x /float(args.n) )
                if cycle==1:  genesPerSample=tmpg # create genesPerSample at first cycle 
                else: genesPerSample=genesPerSample.join(tmpg, on='gene_symbol', how='outer', lsuffix='_genesPerSample', rsuffix='_tmpg').fillna(0)

        #variantsPerGene.to_csv('ciccivar', sep='\t', index=False)    
        genesPerSample['GrandMean']=genesPerSample.sum(axis=1 )/float(args.i)
        genesToDiscard=genesPerSample[genesPerSample['GrandMean']> float(args.gt)]
        #genesPerSample.to_csv('ciccigene', sep='\t', index=False)
        #genesToDiscard['gene_symbol'].to_csv('ciccigenediscard', sep='\t', index=False)

"""
        ###~~~  GREP 
        sample = pd.read_table(args.scsq) 
        sample_filtered=filter(sample, 'genic', args.r, args.pli, args.g, args.cadd)

        listSamples = [line.rstrip('\n') for line in open(args.sl)]
        mapperAFsamples = {}
        for key in sample_filtered.index_x.unique():
            (chrom,pos,alternate) = key.split(":")
            genotypesToConsider, refAllele,  altAlleles = makegenotypelist(args.rvcf, key , listSamples)
            mapperAFsamples.update(Freq_CSQ_REF_ALT (csqAllele, refAllele, altAlleles, missing_data_format, genotypeslist))    
        df_info["AF"] = control_filtered.index_x.map(mapperAF)



        for gid in  listGREP: 

            df= mergeThis(args.f, gid)
    	
            ##~~~  filter in genic regions
            #~1. filter according to criteria 
            df_filtered_g = df[ (df['type']== 'genic') & (df.impact!="MODIFIER") & (df.rare != args.r) & (df.pLI >= args.pli) & ( df.sumGene >= args.g) & (df.CADD >= df.CADD.quantile(args.cadd)) & (df.csqCount >= args.count)]

        #~2. check genes with too many variants 
        variantsPerGene=df_filtered_g[['index_x' , 'gene_symbol']].drop_duplicates().groupby(['gene_symbol']).count()

#####################################################################
        ################### TO ADD, EXCLUDE GENESTODISCARD 
#####################################################################

        df_filtered_g.to_csv(args.o, sep="\t",index=True)

        ##~~~  filter in regulatory regions
        #df_regulatory=df.loc[df['type'] =='regulatory'] 

        ##~~~   filter in intergenic 
"""


if __name__ == "__main__":
        main()
