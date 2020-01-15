# -*- coding: utf-8
import re, sys, argparse, gzip, requests, json 
sys.path.append('/home/enza/ezcngit/grep/filtering/greplib.py')
import greplib as gp


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
    dSOTermFineRank=VepRankingInfo(args.v)
         
    ##### 0b. read weights 
    dWeig={}
    for wline in open (args.w):
        w=wline.rstrip().split() 
        dWeig[w[1]]=int(w[2])
    #print (dWeights)
    
    ##### 1. parse vcf to produce dVcf[mykey]=[myref, myqual, dFormat["GT"]]; mykey is  1:333333:/T (T is the alternate allele) "
    #filemyres=open(args.o, 'w')
    listOfErrors=[]
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
                            dVcf[mykey]=[myref, myqual, dFormat["GT"]]

    ##### 2. load lists of genes          
    lethalFile=open('lethal_candidate.bed','r')
    lethalList=[]
    for i in lethalFile: lethalList.append(i.strip('\n'))
    
    ##### 3. get VEP info 
    dVep={}
    for locusID in dVcf.keys(): 
        #print(locusID)
        dVepValue=getfreqfromVep (locusID)
        #print("ho finito dVep")
        print(dVepValue) 
        if dVepValue: 
            if 'gene_symbol' in dVepValue.keys():
                #print("guardo la lista")
                lethalValue=checkInList(dVepValue, lethalList)
                dVepValue["lethal"]=lethalValue
      
            dVep[locusID]=dVepValue
        else: 
            listOfErrors.append(locusID)
    print(dVep)
    print(dVcf)
    fileToWrite=open(args.e, 'w')
    for i in listOfErrors: fileToWrite.write( i )
 

if __name__ == "__main__":
    main()
