# -*- coding: utf-8
import re, sys, argparse, gzip, json  #requests 
sys.path.append('../libraries')
import greplib as gp
import pandas as pd
import numpy as np
#print("ho caricato i moduli")


#########################################################

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="path to  input  file ",type=str,required=True)
    parser.add_argument("-j", help="path to  json  file ",type=str,required=True)
    parser.add_argument("-g", help="path to gene file list ",required=True)
    #parser.add_argument("-i", help="threshold for SOTerm Impact  ", type=int,required=True)
    parser.add_argument("-r", help="threshold for rare variant definition ", type=float,required=True) 
    parser.add_argument("-v", help="path to table of vep consequences  ",type=str, required= True)   
    parser.add_argument("-o", help="path to output file  ",type=str, required= True)
    #parser.add_argument("-st", help="path to stats file  ",type=str, required= True)
    parser.add_argument("-e", help="path to error file",type=str,required=True)
    parser.add_argument("-c", help="consequence allele count ",type=int,required=False, default=0)
    parser.add_argument("-w", help="path to  weights  file ",type=str,required=True)
    args = parser.parse_args()
    #output = open(args.o,'w')
    #print(args)
        
    ##### 0a. retrieve VEP ranking info   
    dSOTermFineRank=gp.VepRankingInfo(args.v)
         
    ##### 0b. read weights 
    dWeig={}
    for wline in open (args.w):
        w=wline.rstrip().split() 
        dWeig[w[1]]=int(w[2])
    #print (dWeights)
    
    ##### 1. parse vcf to produce dVcf[mykey]=[myref, myqual, dFormat["GT"]]; mykey is  1:333333:/T (T is the alternate allele) "
    
    dVcf={}
    for line in gzip.open(args.f, 'r'):
        decodedLine=line.decode()  ## line.decode() is necessary to read encoded data using gzip in python3
        if not re.match('#', decodedLine):
            linesplit=decodedLine.rstrip().split()
            mychr=linesplit[0]; mypos=linesplit[1]; myref=linesplit[3]; myalt=linesplit[4]; myqual=float(linesplit[5]); altAlleles=myalt.split(",")
            tempformattitle=linesplit[8].split(":")
            tempformatcontent=linesplit[9].split(":")
            dFormat=dict(zip(tempformattitle, tempformatcontent))

            for altAl in altAlleles:
                mykey=mychr.lstrip("chr") + ":" + mypos + ":/" + altAl
                dVcf[mykey]=[myref, myalt, myqual, dFormat["GT"]]
    #print("ho letto il vcf")
   
        
    ##### 2. load genes list
    gene_list = pd.read_csv(args.g,sep="\t")
    #print("ho letto la lista dei geni")

    ##### 3. get VEP info 
    dVep = gp.getInfoFromVepLocally (args.j , args.r )	
    
    """
    ###### 4. create dataframe from dVEP
    '''
                     most_severe_consequence            id csqAllele          gene_id gene_symbol gnomad_nfe gnomad_eas gnomad_asj ....
    5:64548:/A        intergenic_variant           NaN       NaN              NaN         NaN        NaN        NaN        NaN ....
    1:10623:/C     upstream_gene_variant  rs1207502826         C  ENSG00000223972     DDX11L1        NaN        NaN        NaN ....
    1:95068:/A            intron_variant   rs867757274         A  ENSG00000238009  AL627309.1        NaN        NaN        NaN ....
    1:942451:/C         missense_variant     rs6672356         C  ENSG00000187634      SAMD11     0.9999          1          1 ....
    1:1519044:/T          intron_variant   rs147532057         T  ENSG00000197785      ATAD3A        NaN        NaN        NaN ....
    '''

    df = pd.DataFrame(dVep).T
    print (df) 

    #print("ho creato il DataFrame")
    ####### 5 check for common genes and add it!
    common_gene = set(df.gene_id).intersection(set(gene_list.ensID))
    df.loc[:,"score_gene_list"] = df.gene_id.apply(lambda x: gene_list[gene_list.ensID == x].final_score.sum())
    gene_list.loc[:,"value"] = 1
    pivot_tmp = pd.pivot_table(columns="gene_type",index="ensID",data=gene_list,values="value").fillna(0).reset_index()
    df = (df.reset_index().merge(pivot_tmp,how="left",left_on="gene_id",right_on="ensID").drop("ensID",axis=1)).rename({"index":"variant"},axis=1)

    col_freq = ['afr', 'amr', 'gnomad_oth', 'gnomad_fin', 'gnomad', 'gnomad_eas', 'sas','gnomad_afr', 'eur', 'eas', 'gnomad_amr', 'gnomad_asj', 'gnomad_sas','gnomad_nfe']
    common = set(df.columns).intersection(col_freq) 

    def label_race (row):
        if (row[common] < args.r).any(): ## TODO: instead of 0.05 use args.r 
            return 1
        if row[common].isnull().all():
            return 0.5
        return 0

    df.loc[:,"rare"] = df.apply(lambda row: label_race(row), axis=1)

    soScore = pd.Series(dSOTermFineRank,name="soScore").to_frame().reset_index()
    df_last = df.merge(soScore,left_on="most_severe_consequence",right_on="index").drop("index",axis=1)
    #drop columns
    df_last = df_last.set_index("variant").drop(common,axis=1)
    df_last[gene_list.gene_type.unique()] = df_last[gene_list.gene_type.unique()].fillna(0)

    #df_for_stats = df_last[['most_severe_consequence', 'csqCount','csqAllele','gene_id', 'score_gene_list', 'wCellCycle', 'wDDD','wEmbryoDev', 'wEssential', 'wLethal', 'wMisc', 'rare', 'soScore']]

    #final score
    df_last.loc[:,"gpScore"] = (df_last.csqCount.astype(float) * dWeig[wCAC]) + (df_last.rare.astype(float) * dWeig[wRare]) + (df_last.soScore.astype(float) * dWeig[wRank]) + df_last.score_gene_list.astype(float)

    df_last.to_csv(args.o,sep="\t",index=True)
    #df_for_stats.to_csv(args.st,sep="\t",index=True)


if __name__ == "__main__":
    main()
