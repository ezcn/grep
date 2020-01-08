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

def getfreqfromVep (Position):
    freq_dict={}
    server="https://rest.ensembl.org"
    ext = "/vep/human/region/"+ Position +"?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded= r.json()
    if "colocated_variants" in decoded[0]:
        if "id" in decoded[0]["colocated_variants"][0] :
            freq_dict["id"]=decoded[0]["colocated_variants"][0]["id"]
        for var in decoded[0]["colocated_variants"] :
            if "frequencies" in var: freq_dict=var["frequencies"]

    if "most_severe_consequence" in decoded[0]:
        freq_dict["most_severe_consequence"]=decoded[0]["most_severe_consequence"]
        most=freq_dict["most_severe_consequence"]
        if 'transcript_consequences' in decoded[0]: 
            for i in decoded[0]['transcript_consequences']:
                if most in  i['consequence_terms'] :
                    csqAllele=i['variant_allele']
                    freq_dict['csqAllele']=csqAllele
                    freq_dict['gene_id']=i['gene_id']
                    freq_dict['gene_symbol']=i['gene_symbol']
                else:
                    if 'regulatory_feature_consequences' in decoded[0]: 
                        for r  in decoded[0]['regulatory_feature_consequences']:
                            if most in  r['consequence_terms']: 
                                csqAllele=r['variant_allele']
                                freq_dict['csqAllele']=csqAllele
            return freq_dict


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

    ## READ weights 
    dWeig={}
    for wline in open (args.w):
        w=wline.rstrip().split() 
        dWeig[w[1]]=int(w[2])
    #print (dWeights)
        ##  READ VEP consequences rank ########
    """read external file with info on VEP consequences  """
    dRank={"HIGH":4, "LOW": 2, "MODERATE":3, "MODIFIER":1}
    dSOTermRank={}
    lSOTerm=[]  ### list of SOTerm ordered by severity

    countlinesCsq= True
    for csqLine in open(args.v, 'r'):
        if countlinesCsq:
            csqTitle=csqLine.rstrip().split('\t')
            countlinesCsq=False
        else:
            myRowList=csqLine.rstrip().split('\t')
            dCsq= dict(zip(csqTitle, myRowList ))
            dSOTermRank[dCsq['SO term']]=dRank[dCsq['IMPACT']]
            lSOTerm.append(myRowList[0])

    #print (lSOTerm)
    lScores=list(reversed(range(len(lSOTerm)))) 
    #print (lScores) 
    dSOTermFineRank=dict(zip(lSOTerm, map(int, lScores) ))
    #print (dSOTermFineRank)


##########~~~~~~~~~~~~~~  Loop of vcf lines 
    #filemyres=open(args.o, 'w')
    listOfErrors=[]
    ##### 1. parse vcf 
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
        print(locusID)
        dVepValue=getfreqfromVep (locusID)
        print("ho finito dVep")
        print(dVepValue) 
        if dVepValue: 
            if 'gene_symbol' in dVepValue.keys():
                print("guardo la lista")
                lethalValue=checkInList(dVepValue, lethalList)
                dVepValue["lethal"]=lethalValue
      
            dVep[locusID]=dVepValue
        else: 
            listOfErrors.append(locusID)
    print(dVep) 
    fileToWrite=open(args.e, 'w')
    for i in listOfErrors: fileToWrite.write( i )
 

if __name__ == "__main__":
    main()
