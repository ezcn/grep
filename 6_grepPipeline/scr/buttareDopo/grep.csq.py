#!/usr/bin/python3
#import grepLib as gp
import pandas as pd
import numpy as np
import re, sys, argparse, gzip, json, subprocess

###### Functions START ######################## 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def tabix_cadd(key, db ):
    #db = "/lustre/home/enza/CADD/whole_genome_SNVs.tsv.gz"
    (chrom,pos,alternate) = key.split(":")
    cmd = "tabix %s %s:%d-%d" % (db, chrom, int(pos), int(pos))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    result, err = proc.communicate()
    if err: raise IOError("** Error running %s key for %s on %s" % (keyString, db))
    mapper = {}
    cadd_out = result.decode().split("\n")[0:-1]
    for x in cadd_out:
        chrom,pos,ref,alt,score,tmp = x.split("\t")
        mapper[chrom+":"+pos+":/"+alt] = score
    return mapper

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def tabix_fathmm(key, db ):
    #db = "/lustre/home/enza/CADD/whole_genome_SNVs.tsv.gz"
    (chrom,pos,alternate) = key.split(":")
    cmd = "tabix %s %s:%d-%d" % (db, chrom, int(pos), int(pos))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    result, err = proc.communicate()
    if err: raise IOError("** Error running %s key for %s on %s" % (keyString, db))
    mapper = {}
    cadd_out = result.decode().split("\n")[0:-1]
    for x in cadd_out:
        chrom,pos,ref,alt,score = x.split("\t")
        mapper[chrom+":"+pos+":/"+alt] = score
    return mapper

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getInfoFromVepLocally (jsonWithVEPannotations):
    vepInfo={}; vepInfoCommon={}    
    #~~ filename is the json output of VEP runned locally 
    #~~ filename is parsed by line and by alternate allele  
    with open(jsonWithVEPannotations, 'r') as f:
        for line in f:
            info=json.loads(line) #~~ make a dictionary out of a string 
            #~~ derive the key chr:pos:/alt from the "input" element and process each alternate allele 
            locusdata=info['input'].split(); altAlleles=locusdata[4].split(","); mychr=locusdata[0]; mypos=locusdata[1]; mygen=locusdata[9].split(':')[0]
            for altAl in altAlleles:
                mykey=mychr.lstrip("chr") + ":" + mypos + ":/" + altAl
                #print(mykey) 
                most=info["most_severe_consequence"]
                vepInfoCommon[mykey]={}
                vepInfoCommon[mykey]['genotype']=mygen
                vepInfoCommon[mykey]['most_severe_consequence']=most 
                #~~  check if the most sequence is in a transcript
                if 'transcript_consequences' in info:
                    for tc in info['transcript_consequences']:
                        if most in  tc['consequence_terms'] :
                            tcTranscript=tc['transcript_id']
                            vepInfo[(mykey, tcTranscript)]={}
                            tcAllele=tc['variant_allele']
                            if tcAllele ==altAl :
                                vepInfo[(mykey, tcTranscript)]['csqAllele']=tcAllele
                                vepInfo[(mykey, tcTranscript)]['gene_id']=tc['gene_id']
                                vepInfo[(mykey, tcTranscript)]['gene_symbol']= tc['gene_symbol']
                                vepInfo[(mykey, tcTranscript)]['impact']=tc['impact']
                                vepInfo[(mykey, tcTranscript)]['key']=mykey
                                vepInfo[(mykey, tcTranscript)]['element_id']=tc['transcript_id']
                                vepInfo[(mykey, tcTranscript)]['type']='genic'
                #~~ check if the most severe consequence is in a regulatory feature  
                elif 'regulatory_feature_consequences' in info:
                    for rf  in info['regulatory_feature_consequences']:
                        if most in  rf['consequence_terms']:
                            rfRegulatory=rf['regulatory_feature_id']
                            vepInfo[(mykey, rfRegulatory)]={}
                            rfAllele=rf['variant_allele']
                            if rfAllele ==altAl :
                                vepInfo[(mykey, rfRegulatory)]['csqAllele']=rfAllele
                                vepInfo[(mykey, rfRegulatory)]['impact']=rf['impact']
                                vepInfo[(mykey, rfRegulatory)]['key']=mykey
                                vepInfo[(mykey, rfRegulatory)]['element_id']=rf['regulatory_feature_id']
                                vepInfo[(mykey, rfRegulatory)]['type']='regulatory'
                #~~ check if the most severe consequence is in an intergenic region 
                elif 'intergenic_consequences' in info:
                    for ic in info['intergenic_consequences']:
                        if most in  ic['consequence_terms']:
                            vepInfo[(mykey, 'intergenic')]={}
                            icAllele=ic['variant_allele']
                            if icAllele ==altAl:    
                                vepInfo[(mykey, 'intergenic') ]['csqAllele']=icAllele
                                vepInfo[(mykey, 'intergenic') ]['impact']=ic['impact']
                                vepInfo[(mykey, 'intergenic') ]['key']=mykey
                                vepInfo[(mykey, 'intergenic') ]['element_id']='intergenic'
                                vepInfo[(mykey, 'intergenic') ]['type']='intergenic'
                #~~ retrive info on rsID, starting position, frequencies for the csqAllele found in the previous code  
                if "colocated_variants" in info:
                    infoCV = info["colocated_variants"][0]
                    #~~ adding info for starting position, this is needed for merging with CADD score.
                    vepInfoCommon[mykey]["start"] = infoCV["start"]

                    #~~ check if rsid is present 
                    if "id" in infoCV: vepInfoCommon[mykey]["id"] = infoCV["id"]
                    
                    #~~ retrieve allelle frequencies 
                    if 'frequencies' in infoCV:
                        vepInfoCommon[mykey]["frequencies"]= infoCV["frequencies"]
                else:
                    infoCV={}                          
    return vepInfo, vepInfoCommon
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def csqAlleleCount (csqAllele, altAllele, altAlleleCount, ploidy): 
    #csqAllele= cons allele  from vep  
    #altAlleleCount= integer, alternate allele count 
    # ploidy = number of allels at one locus 
    
    if csqAllele == altAllele: mycsqAlleleCount = altAlleleCount      ## Consequence allele is the alternate allele 
    else: mycsqAlleleCount = ploidy-altAlleleCount    ## Consequence allele is the reference allele;
    
    return mycsqAlleleCount
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def checkFreq (listFreq, threshold): 
    # check if a variant is rare:  none of the populations in listfreq has frequency greater than threshold
    rare=True
    if len(listFreq) >0: 
        for freq in listFreq:
            if freq > threshold:
                rare=False
                break
    else: rare="NOB" #never observed 
    return(rare)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def VepRankingInfo (vepinfofile): 
    """read external file with info on VEP consequences  """
    dRank={"HIGH":4, "LOW": 2, "MODERATE":3, "MODIFIER":1}
    dSOTermRank={}
    lSOTerm=[]  ### list of SOTerm ordered by severity

    countlinesCsq= True
    for csqLine in open(vepinfofile, 'r'):
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
    return dSOTermFineRank

###### Functions END ###################################


def main():
    parser = argparse.ArgumentParser()
    #~~~ input files 
    parser.add_argument("-j", help="path to  json  file ",type=str,required=True)
    parser.add_argument("-g", help="path to gene file list ",required=True)
    
    parser.add_argument("-pli", help="path to table of pLI score  ",type=str, required= True)   
    parser.add_argument('-cadd', help="path to CADD file",type=str, required= True)   
    parser.add_argument('-fathmm', help="path to fathmm file",type=str, required= True)   
    parser.add_argument('-fathmmcoding', help="path to fathmm coding file",type=str, required= True)   

    parser.add_argument("-v", help="path to table of vep consequences  ",type=str, required= True)

    parser.add_argument("-c", help="consequence allele count ",type=int,required=False, default=0)
    
    parser.add_argument("-r", help="threshold for rare variant definition ", type=float,required=True) 
    
    parser.add_argument("-o", help="path to output file  ",type=str, required= True)
    #parser.add_argument("-st", help="path to stats file  ",type=str, required= True)
    parser.add_argument("-e", help="path to error file",type=str,required=True)

    args = parser.parse_args()
        
    
    ##### 1. parse VEP info from local json file  
    dV = getInfoFromVepLocally (args.j )   #yelds two dictionaries 
    dVepTrans=dV[0]; dVepCommon=dV[1]

    #~~ at each  locus estimate  the number of allele with consequences 
    # TODO: single step, apply function in pandas    
    for telem in dVepTrans:
        if 'csqAllele' in dVepTrans[telem]:     
            vcfkey=telem[0]; csqAllele=dVepTrans[telem]['csqAllele']; altAllele=vcfkey.split('/')[-1]; genotype=dVepCommon[vcfkey]['genotype']
            altAlleleCount=2-genotype.count('0')
            csqCount= csqAlleleCount (csqAllele, altAllele, altAlleleCount, 2)
            dVepTrans[telem]['csqCount']=csqCount  #1 het #2 homo

    #~~check if variant is rare 
    for celem in dVepCommon:
        if 'frequencies' in dVepCommon[celem]: 
            for cAll in dVepCommon[celem]['frequencies']: 
                listFreq=dVepCommon[celem]['frequencies'][cAll].values()
                rare=checkFreq (listFreq, args.r )
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
    dSOTermFineRank=VepRankingInfo(args.v)
    soScore = pd.Series(dSOTermFineRank,name="soScore").to_frame().reset_index()
    df_final = df_last.merge(soScore,left_on="most_severe_consequence",right_on="index").set_index("index_x").drop("index_y",axis=1)
    #df_final.to_csv(args.o,sep="\t",index=True)

    #### 4. pLI 
    ### load pli table and add pli to df
    pLI_score = pd.read_csv(args.pli,sep="\t")
    pliScore=pLI_score[["transcript", "pLI"]]
    df_info = df_final.reset_index().merge(pliScore,left_on="element_id",right_on="transcript",how="left").drop(["transcript", "frequencies", "genotype"] , axis=1)

    #### 5. CADD
    mappercadd = {}
    for key in df_info.index_x.unique():
        mappercadd.update(tabix_cadd(key, args.cadd))

    df_info["CADD"] = df_info.index_x.map(mappercadd)

    #### 6. fathmm-xf non-coding 
    mapperfathmm = {}
    for key in df_info.index_x.unique():
        mapperfathmm.update(tabix_fathmm(key, args.fathmm))

    df_info["fathmm-xf"] = df_info.index_x.map(mapperfathmm)

    #### 6. fathmm-xf coding 
    mapperfathmmCoding = {}
    for key in df_info.index_x.unique():
        mapperfathmmCoding.update(tabix_fathmm(key, args.fathmmcoding))

    df_info["fathmm-xfCoding"] = df_info.index_x.map(mapperfathmmCoding)


    #### 7. print output tsv 
    df_info.to_csv(args.o,sep="\t",index=False)
  

if __name__ == "__main__":
    main()
