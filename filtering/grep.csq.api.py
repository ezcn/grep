# -*- coding: utf-8
import re, sys, argparse, gzip, requests, json 
sys.path.append('/home/enza/ezcngit/grep/filtering/greplib.py')
import greplib as gp
import pandas as pd
import numpy as np



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

    ##### 2. load genes list
    gene_list = pd.read_csv(args.m,sep="\t")
    
    ##### 3. get VEP info 
    listOfErrors=[]
    dVep={}
    for locusID in dVcf.keys(): 
        #print(locusID)
        dVepValue=gp.getInfoFromVep (locusID)
        #print("ho finito dVep")
        #prnt(dVepValue) 
        if dVepValue: 
            dVep[locusID]=dVepValue
            dVep["csqCount"]= gp.CountCSQ_REF_ALT(dVep["csqAllele"], dVcf[locusID][0], dVcf[locusID][1], [dVcf[locusID][3]]) [0]
            
        else: 
            listOfErrors.append(locusID)

    filemyres=open(args.o, 'w')
    for vv in dVep: vvstring= vv + ' # ' + str(dVep[vv]) + '\n'; filemyres.write(vvstring)
    for vc in dVcf: vcstring= vc + ' # ' + str(dVcf[vc]) + '\n' ; filemyres.write(vcstring)

    fileToWrite=open(args.e, 'w')
    for i in listOfErrors: fileToWrite.write( i )
    
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

    ####### 5 check for common genes and add it!
    common_gene = set(df.gene_id).intersection(set(gene_list.ensID))
    df.loc[:,"score_gene_list"] = df.gene_id.apply(lambda x: gene_list[gene_list.ensID == x].final_score.sum())




if __name__ == "__main__":
    main()
