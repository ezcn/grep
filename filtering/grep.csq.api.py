# -*- coding: utf-8
import re, sys, argparse, gzip, requests, json 
sys.path.append('/lustrehome/gianluca/scripts/greplib.py')
import greplib as gp
import pandas as pd
import numpy as np
print("ho caricato i moduli")


########################################################

def checkInList(gene, listOfGenes):
        #gene=(mydict['gene_symbol'])
        IsIn=False
        if gene in listOfGenes:
            IsIn=True
        return IsIn
#########################################################

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="path to  input  file ",type=str,required=True)
    parser.add_argument("-g", help="path to gene file list ",required=True)
    #parser.add_argument("-i", help="threshold for SOTerm Impact  ", type=int,required=True)
    parser.add_argument("-r", help="threshold for rare variant definition ", type=float,required=True) 
    parser.add_argument("-v", help="path to table of vep consequences  ",type=str, required= True)   
    parser.add_argument("-o", help="path to output file  ",type=str, required= True)
    parser.add_argument("-e", help="path to error file",type=str,required=True)
    parser.add_argument("-c", help="consequence allele count ",type=int,required=False, default=0)
    parser.add_argument("-w", help="path to  weight  file ",type=str,required=True)
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
    print("ho letto il vcf")
    ##### 2. load genes list
    gene_list = pd.read_csv(args.g,sep="\t")
    print("ho letto la lista dei geni")
    ##### 3. get VEP info 
    listOfErrors=[]
    dVep={}
    for locusID in dVcf.keys(): 
        #print(locusID)
        dVepValue=gp.getInfoFromVep (locusID)
        #print("ho finito dVep")
        #prnt(dVepValue) 
        if type(dVepValue) is not str: 
            dVep[locusID]=dVepValue
            if "csqAllele" in dVep[locusID]:
                dVep[locusID]["csqCount"]= gp.CountCSQ_REF_ALT(dVep[locusID]["csqAllele"], dVcf[locusID][0], dVcf[locusID][1], [dVcf[locusID][3]]) [0]
            else:
                dVep[locusID]["csqCount"] = np.nan
        
        else: 
            listOfErrors.append(locusID)

    #filemyres=open(args.o, 'w')
    #for vv in dVep: vvstring= vv + ' # ' + str(dVep[vv]) + '\n'; filemyres.write(vvstring)
    #for vc in dVcf: vcstring= vc + ' # ' + str(dVcf[vc]) + '\n' ; filemyres.write(vcstring)
    print("ho preso le info da VEP")
    fileToWrite=open(args.e, 'w')
    for i in listOfErrors: fileToWrite.write( i )
    print("ho scritto il file di errore")
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
    print("ho creato il DataFrame")
    ####### 5 check for common genes and add it!
    common_gene = set(df.gene_id).intersection(set(gene_list.ensID))
    df.loc[:,"score_gene_list"] = df.gene_id.apply(lambda x: gene_list[gene_list.ensID == x].final_score.sum())

    pop = ['afr', 'amr', 'gnomad_oth', 'gnomad_fin', 'gnomad', 'gnomad_eas','sas','gnomad_afr', 'eur', 'eas', 'gnomad_amr', 'gnomad_asj', 'gnomad_sas','gnomad_nfe']

    def label_race (row):
        if (row[pop] < 0.05).any():
            return 1
        if row[pop].isnull().all():
            return 0.5
        return 0

    df.loc[:,"rare"] = df.apply(lambda row: label_race(row), axis=1)

    soScore = pd.Series(dSOTermFineRank,name="soScore").to_frame().reset_index()
    df_last = df.reset_index().merge(soScore,left_on="most_severe_consequence",right_on="index").rename({"index_x" : "variant"},axis=1).drop("index_y",axis=1)
    #drop columns
    df_last.drop(pop,axis=1,inplace=True)

    wCAC = 1
    wRare = 1
    wRank = 1
    #final score
    df_last.loc[:,"gpScore"] = (df_last.csqCount.astype(float) * wCAC) + (df_last.rare.astype(float) * wRare) + (df_last.soScore.astype(float) * wRank) + df_last.score_gene_list.astype(float)

    df_last.to_csv(arg.o,sep="\t",index=False)

if __name__ == "__main__":
    main()
