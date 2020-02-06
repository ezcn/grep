import re, sys, argparse, gzip, requests, json
import pandas as pd
import numpy as np

dVcf={}
for line in gzip.open("/lustrehome/giuliana/mity/AS006.mity.vcf.gz", 'r'):
	decodedLine=line.decode()  ## line.decode() is necessary to read encoded data using gzip in python3
	if not re.match('#', decodedLine):
		linesplit=decodedLine.rstrip().split()
		mychr=linesplit[0]; mypos=linesplit[1]; myref=linesplit[3]; myalt=linesplit[4]; myqual=float(linesplit[5]); altAlleles=myalt.split(",")
		tempformattitle=linesplit[8].split(":")
		tempformatcontent=linesplit[9].split(":")
		dFormat=dict(zip(tempformattitle, tempformatcontent))
		for altAl in altAlleles:
			mykey="chr" + mychr.lstrip("chr") + ":" + mypos + ":/" + altAl
			dVcf[mykey]=[myref, myalt, myqual, dFormat["GT"]]




def getInfoFromVep (Position):
	mitidic={}
	server="https://rest.ensembl.org"
	ext = "/vep/human/region/"+ Position +"?"
	r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
	if not r.ok:
		return Position
		#r.raise_for_status()
		#sys.exit()
	decoded= r.json()
	info = decoded[0]
	if "colocated_variants" in info:
		id_search = info["colocated_variants"][0]
		if "id" in id_search :
			mitidic["id"] = id_search["id"]
	if "most_severe_consequence" in info:
		mitidic["most_severe_consequence"]=info["most_severe_consequence"]
		most=mitidic["most_severe_consequence"]
		if 'transcript_consequences' in info: 
			for i in info['transcript_consequences']:
				if most in  i['consequence_terms'] :
					csqAllele=i['variant_allele']
					mitidic['csqAllele']=csqAllele
					mitidic['gene_id']=i['gene_id']
					mitidic['gene_symbol']=i['gene_symbol']
					mitidic['impact']=i['impact']
				else:
					if 'regulatory_feature_consequences' in info: 
						for r  in info['regulatory_feature_consequences']:
							if most in  r['consequence_terms']: 
								csqAllele=r['variant_allele']
								mitidic['csqAllele']=csqAllel
								mitidic['impact']=r['impact'] 
	return mitidic


dVep={}
for Position in dVcf.keys(): 
	dVepValue=getInfoFromVep (Position)
	if dVepValue: 
		dVep[Position]=dVepValue

df = pd.DataFrame(dVep).T

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
    lScores=list(reversed(range(len(lSOTerm)))) 
    #print (lScores) 
    dSOTermFineRank=dict(zip(lSOTerm, map(int, lScores) ))
    #print (dSOTermFineRank)
    return dSOTermFineRank
    
dSOTermFineRank=VepRankingInfo("/lustrehome/giuliana/github/grep/filtering/csqimpact.tsv")
SoScore = pd.Series(dSOTermFineRank,name="soScore").to_frame().reset_index()
df=df.reset_index().merge(SoScore,left_on="most_severe_consequence",right_on="index").set_index("index_x").drop("index_y",axis=1)

df.to_csv("/lustrehome/giuliana/mity/VEP_mt/test_impact",sep="\t",index=True)



