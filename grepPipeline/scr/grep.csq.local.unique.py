#!/usr/bin/python3
import greplib as gp
import pandas as pd
import numpy as np
import re, sys, argparse, gzip, json, subprocess
import dask.dataframe as dd
from ast import literal_eval
import ast
import decimal



def main():
    parser = argparse.ArgumentParser()
    #~~~ input files 
    parser.add_argument("-j", help="path to  json  file ",type=str,required=True)
    parser.add_argument("-g", help="path to gene file list ",required=True)
    
    parser.add_argument("-pli", help="path to table of pLI score  ",type=str, required= True)   
    parser.add_argument('-cadd', help="path to CADD file",type=str, required= True)   
    parser.add_argument('-fathmm', help="path to fathmm file",type=str, required= True)   

    parser.add_argument("-v", help="path to table of vep consequences  ",type=str, required= True)

    parser.add_argument("-c", help="consequence allele count ",type=int,required=False, default=0)
    
    parser.add_argument("-r", help="threshold for rare variant definition ", type=float,required=True) 
    
    parser.add_argument("-o", help="path to output file  ",type=str, required= True)
    #parser.add_argument("-st", help="path to stats file  ",type=str, required= True)
    parser.add_argument("-e", help="path to error file",type=str,required=True)

    args = parser.parse_args()
        
    
    ##### 1. parse VEP info from local json file  
    dV = gp.getInfoFromVepLocally (args.j )   #yelds two dictionaries 
    dVepTrans=dV[0]; dVepCommon=dV[1]

    #~~ at each  locus estimate  the number of allele with consequences 
    # TODO: single step, apply function in pandas    
    for telem in dVepTrans:
        if 'csqAllele' in dVepTrans[telem]:     
            vcfkey=telem[0]; csqAllele=dVepTrans[telem]['csqAllele']; altAllele=vcfkey.split('/')[-1]; genotype=dVepCommon[vcfkey]['genotype']
            altAlleleCount=2-genotype.count('0')
            csqCount= gp.csqAlleleCount (csqAllele, altAllele, altAlleleCount, 2)
            dVepTrans[telem]['csqCount']=csqCount  #1 het #2 homo

    #~~check if variant is rare 
    for celem in dVepCommon:
        if 'frequencies' in dVepCommon[celem]: 
            for cAll in dVepCommon[celem]['frequencies']: 
                listFreq=dVepCommon[celem]['frequencies'][cAll].values()
                rare=gp.checkFreq (listFreq, args.r )
                dVepCommon[celem]['rare']=rare
                dVepCommon[celem]['commonCSQall']=cAll
        else:
                dVepCommon[celem]['rare']="NOB"
 
    #~~create dataframe from dV
    dfTrans = pd.DataFrame(dVepTrans).T 
    dfCommon= pd.DataFrame(dVepCommon).T
    #dfTrans.to_csv("cicci.trans.tsv", sep="\t")#, index=True )
    #dfCommon.to_csv("cicci.common.tsv", sep="\t")#, index=True )
    
    df=dfTrans.set_index('key').join(dfCommon, how="left"  )
    #df.to_csv("/lustrehome/silvia/cicci.ttcos.tsv", sep="\t")
    
    #### 2. info form gene lists
    gene_list = pd.read_csv(args.g,sep="\t")
    df_last=df.reset_index().merge(gene_list, left_on="gene_id", right_on="ensID", how="left").drop("ensID", axis=1).fillna(0)
 
    #### 3. SOScore 
    dSOTermFineRank=gp.VepRankingInfo(args.v)
    soScore = pd.Series(dSOTermFineRank,name="soScore").to_frame().reset_index()
    df_final = df_last.merge(soScore,left_on="most_severe_consequence",right_on="index").set_index("index_x").drop("index_y",axis=1)
    #df_final.to_csv(args.o,sep="\t",index=True)

    #### 4. pLI 
    ### load pli table and add pli to df
    pLI_score = pd.read_csv(args.pli,sep="\t")
    pliScore=pLI_score[["transcript", "pLI"]]
    df_info = df_final.reset_index().merge(pliScore,left_on="element_id",right_on="transcript",how="left").drop(["transcript", "frequencies", "genotype"] , axis=1)

    #### 5. CADD
    mapper = {}
    for key in df_info.index_x.unique():
        mapper.update(tabix_cadd(key, args.cadd))

    df_info["CADD"] = df_info.index_x.map(mapper)

    df_info.to_csv(args.o,sep="\t",index=False)

if __name__ == "__main__":
    main()
