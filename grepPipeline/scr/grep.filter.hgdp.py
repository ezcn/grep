#!/usr/bin/python3
import pandas as pd
import glob, argparse,  subprocess, random 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def makeMapperAF (keylist, vcf, listOfSampleID): 
    mapperAF={}
    for key in keylist:
        alternate=key.split(":")[2].lstrip('/')
        genotypesToConsider, refAllele,  altAlleles = takeInfoFromVcf(vcf, key, listOfSampleID)
        alternateAF = Freq_CSQ_REF_ALT (alternate, refAllele, altAlleles, '.', genotypesToConsider)     
        if alternateAF: mapperAF[key] = alternateAF[0]  
    return mapperAF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def takeInfoFromVcf(vcf, key , samplesToConsider):
    (chrom,pos,tmpalternate) = key.split(":")
    alternate=tmpalternate.lstrip('/')

    cmd_head = "tabix -H %s %s:%d-%d | tail -1 " % (vcf, chrom, int(pos), int(pos))
    proc_head = subprocess.Popen(cmd_head, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    result_head, err = proc_head.communicate()
    if err: raise IOError("** Error running %s key for %s on %s" % (keyString, db))
    vcfhead=result_head.decode().rstrip().split('\t')

    cmd = "tabix  %s %s:%d-%d " % (vcf, chrom, int(pos), int(pos))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    result, err = proc.communicate()
    if err: raise IOError("** Error running %s key for %s on %s" % (keyString, db))
    vcfout =[x.rstrip().split('\t') for  x in result.decode().rstrip().split('\n')]
   
    column2retain=[vcfhead.index(ind) for ind in samplesToConsider]
    genotypesToConsider=[locus[indx].split(":")[0] for indx in column2retain for locus in vcfout if locus[4]==alternate]
    #reference=[locus[3] for locus in vcfout if locus[4]==alternate.lstrip('/')][0]
    reference=''
    for locus in vcfout: 
        if alternate  in locus[4].split(','): reference=locus[3]
    #reference=[locus[3] for locus in vcfout if alternate.lstrip('/')  in locus[4].split(',')][0]
    #dtmp=dict(zip(vcfout[0], vcfout[1]))
    #genotypesToConsider=[dtmp[i].split(':')[0] for i in samplesToConsider]
    return  genotypesToConsider , reference, alternate #vcfhead, vcfout #dtmp,
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def Freq_CSQ_REF_ALT (csqAllele, refAllele, altAlleles, missing_data_format, genotypeslist):
    """csqAllele = type: string, consequence allele  from vep  
    refAllele = type: string, reference allele from variant calling
    altAlleles = type: string, comma-separated list of alternate allele from variant calling
    nbAploidSamples : calculated 
    GTfields = type: list, list of genotypes ["0/0", "0/1", ".", "./."]"""
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
       # freqREF="{0:4.2f}".format(dAllele[refAllele]/nbHaploidSamples)
        freqAlt=[]
        for i in listAlt: 
            freqAltTemp="{0:4.2f}".format(dAllele[i]/nbHaploidSamples)
            freqAlt.append(freqAltTemp)
        #~~~ CSQ 
        if csqAllele in dAllele:
            csqAllCount=dAllele[csqAllele]
            freqCsq="{0:4.2f}".format(csqAllCount/nbHaploidSamples) 
        else: freqCsq='NA'
        myres= [freqCsq, csqAllCount]
    return myres

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def filter (df, genometype, thresold_rare, threshold_pli, threshold_sumgene, threshold_cadd): 
    df_filtered = df[ (df['type']== genometype) & (df.impact!="MODIFIER") & (df.rare != thresold_rare) & (df.pLI >= threshold_pli) & ( df.sumGene >= threshold_sumgene) & (df.CADD >= df.CADD.quantile(threshold_cadd))] #& (df.csqCount >= args.count)]
    return df_filtered

#FUNCTIONS END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def main():
        parser = argparse.ArgumentParser()
        #HGDP arguments
        parser.add_argument("-ccsq", help="path to control csq  ",type=str,required=True)
        parser.add_argument("-cvcf", help="path to reference vcf  ",type=str,required=True)
        parser.add_argument("-cl", help="list of reference id ",type=str,required=True)
        parser.add_argument("-n", help="number of individual to sample ",type=int,required=True)
        parser.add_argument("-i", help="number of iterations ",type=int,required=True) 
        parser.add_argument("-gt", help="threshold for excluding genes ",type=float,required=True) 

        #GREP arguments
        parser.add_argument("-scsq", help="path to control csq  ",type=str,required=True)
        parser.add_argument("-svcf", help="path to reference vcf  ",type=str,required=True)
        parser.add_argument("-sl", help="list of reference id ",type=str,required=True)

        #FILTERING arguments    
        parser.add_argument("-r", help=" False, True ", type=str,required=True)
        parser.add_argument("-pli", help="threshold for  pLI score  ",type=float, required= True)
        parser.add_argument("-g", help=" number of gene lists  ",type=float , required= True)
        parser.add_argument("-cadd", help=" treshold for CADD score  ",type=float , required= True)
        parser.add_argument("-ac", help=" allele count >= of   ",type=int , required= True)
        #parser.add_argument("-count", help=" treshold for csq allele count  ",type=float , required= True)

        parser.add_argument("-o", help="path to output file  ",type=str, required= True)
        parser.add_argument("-e", help="path to error file",type=str,required=True)
        args = parser.parse_args()
        #sys.stdout=open(args.o, 'w')   

        #~~~  HGDP: Random sampling args.i times of args.n individual to figure out genes that show up on average args.gt% times over args.i  iterations. Annotate genesToDiscard in a list and exclude these genes from results 
        ##~~ read and filter annotated info on variable loci in controls  
        control = pd.read_table(args.ccsq) 
        control.loc[:,"sumGene"] = control["EmbryoDev"]+ control["DDD"]+ control["Lethal"]+ control["Essential"]+ control["Misc"]
        #control.to_csv('control')
        control_filtered=filter(control, 'genic', args.r, args.pli, args.g, args.cadd)
        #control_filtered.to_csv('control_filtered_pre')

       
        ##~~ repeat args.i times on args.n samples from controls  
        listControl = [line.rstrip('\n') for line in open(args.cl)]
        cycle=0
        while cycle < args.i: 
                cycle+=1

                ##~~ choose a random sample form controls of size args.n  
                randomSample=random.sample(listControl, args.n)
                #`print (randomSample)

                ##~~ integrate with allele freqency calculated in sample to consider
                keylist=control_filtered.index_x.unique()
                mapperAF=makeMapperAF (keylist, args.cvcf, randomSample)
                control_filtered["af"] = control_filtered.index_x.map(mapperAF)
                #control_filtered.to_csv('control_filtered')
                
                ##~~ integrate with samples ID and csqAlele count 
                mapperAC = {}; mapperSS={}
                control_filtered_allsamples=pd.DataFrame()    

                for ss in randomSample:  
                    tmpdf=control_filtered
                    for key in tmpdf.index_x.unique():
                        alternate=key.split(':')[2].lstrip('/')                 
                        ###~~ csq allele counts for sample ss at locus key 
                        genotypesToConsider, refAllele,  altAlleles = takeInfoFromVcf(args.cvcf, key, [ss])
                        alternateAC = Freq_CSQ_REF_ALT (alternate, refAllele, altAlleles, '.', genotypesToConsider)     
                        if alternateAC: 
                            if alternateAC[1] >= args.ac: 
                                mapperAC[key] = alternateAC[1]   
                                mapperSS[key] = ss
                        tmpdf['ac']= tmpdf.index_x.map(mapperAC)
                        tmpdf['sample'] = tmpdf.index_x.map(mapperSS)
                        ####### tmpdf remove rows with no AC 
                        df = tmpdf[tmpdf['ac'].notna()]
                        df = tmpdf[tmpdf['sample'].notna()]
                        #, tmpdf['sample'].notna()]
                        control_filtered_allsamples=pd.concat([control_filtered_allsamples, df]) 

                #control_filtered_allsamples.to_csv('control_filtered_allsamples', index=False)                
                ##~~ subset for loci with AF>0 in sampleToConsider 
                ##~~ count each gene only once per sample -> gene symbol, in how many samples it is found 
                tmpg=control_filtered_allsamples[[ 'gene_symbol', 'sample']].drop_duplicates().groupby('gene_symbol').count().transform(lambda x: x /float(args.n) )
                #tmpg.to_csv('tmpg%s' % (cycle) , sep='\t') #, index=False)

                if cycle==1:  genesPerSample=tmpg # create genesPerSample at first cycle 
                else: genesPerSample=genesPerSample.join(tmpg, on='gene_symbol', how='outer', lsuffix='_genesPerSample', rsuffix='_tmpg').fillna(0)    
        ##~~ make GrandMean over args.i and discard genes that on average shows up in args.gt individuals over args.i iterations 
        genesPerSample['GrandMean']=genesPerSample.sum(axis=1 )/float(args.i)
        genesToDiscard=genesPerSample[genesPerSample['GrandMean']> float(args.gt)]
        #genesPerSample.to_csv('ciccigene', sep='\t') 
        genesToDiscard.to_csv('ciccigenediscard', sep='\t') 


        ###~~~  GREP
        sample = pd.read_table(args.scsq) 
        sample.loc[:,"sumGene"] = sample["EmbryoDev"]+ sample["DDD"]+ sample["Lethal"]+ sample["Essential"]+ sample["Misc"]
        #sample.to_csv('sample')

        ###~~~ GENIC REGIONS 
        sample_filtered=filter(sample, 'genic', args.r, args.pli, args.g, args.cadd)
        #sample_filtered.to_csv('sample_filtered_pre')


        listSamples = [line.rstrip('\n') for line in open(args.sl)]

        ##~~ integrate with allele freqency calculated in sample to consider
        skeylist=sample_filtered.index_x.unique()
        smapperAF=makeMapperAF (skeylist, args.svcf, listSamples)
        sample_filtered["af"] = sample_filtered.index_x.map(mapperAF)
        #sample_filtered.to_csv('sample_filtered')

        ##~~ integrate with samples ID and csqAlele count 
        smapperAC = {}; smapperSS={}
        sample_filtered_allsamples=pd.DataFrame()    

        for ll in listSamples:   
            stmpdf=sample_filtered
            for skey in stmpdf.index_x.unique():
                        salternate=skey.split(':')[2].lstrip('/')                 
                        ###~~ csq allele counts for sample ll at locus skey 
                        sgenotypesToConsider, srefAllele,  saltAlleles = takeInfoFromVcf(args.svcf, skey, [ll])
                        salternateAC = Freq_CSQ_REF_ALT (salternate, srefAllele, saltAlleles, '.', sgenotypesToConsider)     
                        if salternateAC: 
                            if salternateAC[1] >= args.ac: 
                                smapperAC[key] = salternateAC[1]   
                                smapperSS[key] = ll
                        stmpdf['ac']= stmpdf.index_x.map(smapperAC)
                        stmpdf['sample'] = stmpdf.index_x.map(smapperSS)
                        ####### tmpdf remove rows with no AC 
                        sdf = stmpdf[stmpdf['ac'].notna()]
                        sdf = stmpdf[stmpdf['sample'].notna()]
                        #, tmpdf['sample'].notna()]
                        sample_filtered_allsamples=pd.concat([sample_filtered_allsamples, sdf]) 

        sample_filtered_allsamples.to_csv('sample_filtered_allsamples', index=False)             

        ##~~ remove genes that shows up in the control args.gt of times over args.i iterations  
        sample_filtered_allsamples.drop(sample_filtered_allsamples[sample_filtered_allsamples.gene_symbol.isin(genesToDiscard.gene_symbol) ].index, inplace=True
            )
        
        ####  TO DO ##########################################
        ##~~ remove genes with too many variants 
        #variantsPerGene=sample_filtered_allsamples[['index_x' , 'gene_symbol']].drop_duplicates().groupby(['gene_symbol']).count()
        ####  rimuovere varianti in variantsPerGene Possiamo farlo solo dopo che c'e un file di prova csq decente  


        ##~~~ print very final dataframe 
        #df_filtered_g.to_csv(args.o, sep="\t",index=True)

        ##~~~  filter in regulatory regions
        ##~~~   filter in intergenic 



if __name__ == "__main__":
        main()
