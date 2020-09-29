#!/usr/bin/python3
#import requests, json 
import json 

def getInfoFromVepLocally (jsonWithVEPannotations, thresholdRareAllele ):  
	vepInfo={}
	#~~ filename is the json output of VEP runned locally 
	#~~ filename is parsed by line and by alternate allele  
	with open(jsonWithVEPannotations, 'r') as f:  
		for line in f:
			info=json.loads(line) #~~ make a dictionary out of a string 
			#print(line_d['input']) 
			#~~ derive the key chr:pos:/alt from the "input" element and process each alternate allele 
			locusdata=info['input'].split(); altAlleles=locusdata[4].split(","); mychr=locusdata[0]; mypos=locusdata[1]
			for altAl in altAlleles:
				mykey=mychr.lstrip("chr") + ":" + mypos + ":/" + altAl
				vepInfo[mykey]={}      
				#print ('#################' ) 
				#print ( mykey )
				#~~ check if the most severe consequence is in the line else retrun an empty vepInfo[mykey] 
				if "most_severe_consequence" in info:
					vepInfo[mykey]["most_severe_consequence"]=info["most_severe_consequence"]
					most=vepInfo[mykey]["most_severe_consequence"]
					#~~  check if the most sequence is in a transcript
					if 'transcript_consequences' in info:
						ReportedGeneId=[] ; ReportedGeneSymbol=[]
						for tc in info['transcript_consequences']:
							if most in  tc['consequence_terms'] :
								csqAllele=tc['variant_allele']
								if csqAllele ==altAl : 
									#print (most, csqAllele, tc["transcript_id"], tc['gene_id'])
									vepInfo[mykey]['csqAllele']=csqAllele
									if not tc['gene_id'] in ReportedGeneId: ReportedGeneId.append(tc['gene_id'])
									if not tc['gene_symbol'] in ReportedGeneSymbol:  ReportedGeneSymbol.append(tc['gene_symbol']) 
						vepInfo[mykey]['gene_id']=ReportedGeneId
						vepInfo[mykey]['gene_symbol']=ReportedGeneSymbol

                                        #~~ check if the most severe consequence is in a regulatory feature  
					elif 'regulatory_feature_consequences' in info: 
						for rf  in info['regulatory_feature_consequences']:
							if most in  rf['consequence_terms']:
								csqAllele=rf['variant_allele']
								if csqAllele ==altAl :
									vepInfo[mykey]['csqAllele']=csqAllele
						vepInfo[mykey]['gene_id']=[]
						vepInfo[mykey]['gene_symbol']=[]

                                        #~~ check if the most severe consequence is in an intergenic region 
					elif 'intergenic_consequences' in info: 
						for ic in info['intergenic_consequences']: 
							if most in  ic['consequence_terms']: 
								csqAllele=ic['variant_allele']
								vepInfo[mykey]['csqAllele']=csqAllele	 
						vepInfo[mykey]['gene_id']=[]
						vepInfo[mykey]['gene_symbol']=[]

				#print(vepInfo[mykey]) 
				#~~ retrive info on rsID, starting position, frequencies for the csqAllele found in the previous code  
				if "colocated_variants" in info:
					id_search = info["colocated_variants"][0]
					if "id" in id_search:
						vepInfo[mykey]["id"] = id_search["id"]
					#~~ adding info for starting position, this is needed for merging with CADD score.
					if "start" in id_search: 
						vepInfo[mykey]["starting_position"] = id_search["start"]   
					#~~ check if variant is rare according to threshold given as argument and report True or False in vepInfo[mykey]
					if 'frequencies' in id_search:
						frqInfoDic=id_search["frequencies" ]
						for i in frqInfoDic: 
							if 'csqAllele' in vepInfo[mykey]: 
								if vepInfo[mykey]['csqAllele'] in frqInfoDic: 
									to_parse = list(frqInfoDic[vepInfo[mykey]['csqAllele']].values())
									#print(to_parse)
									vepInfo[mykey]['rare']= checkFreq (to_parse, thresholdRareAllele )
                                                        	#else: vepInfo[mykey]['rare']='na'
		

	#print (json.dumps(info, indent=4 ) )  
	#print (json.dumps(vepInfo, indent=4) ) 
	return vepInfo






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def VepSOTermInfo (vepinfofile): 
    """read external file with info on VEP consequences  """
    lSOTerm=[]  ### list of SOTerm 
    
    countlinesCsq= True
    for csqLine in open(vepinfofile, 'r'):
        if countlinesCsq:
            csqTitle=csqLine.rstrip().split('\t')
            countlinesCsq=False
        else:
            myRowList=csqLine.rstrip().split('\t')
            lSOTerm.append(myRowList[0])

    return  lSOTerm

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getInfoFromVep (Position):
    """Retrieves information from Variant Effect Predictor API
    dependencies: requests, json 
    Position =  1:333333:/T (T is the alternate allele)   """
    freq_dict={}
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
        if "id" in id_search:
            freq_dict["id"] = id_search["id"]
        # adding info for starting position, this is needed for merging with CADD score.
        #if "start" in id_search: 
        #	freq_dict["starting_position"] = id_search["start"]   
        if 'frequencies' in id_search:
            to_parse = list(id_search["frequencies"].values())[0]
            for var in to_parse.keys():
                freq_dict[var]=to_parse[var]

    if "most_severe_consequence" in info:
        freq_dict["most_severe_consequence"]=info["most_severe_consequence"]
        most=freq_dict["most_severe_consequence"]
        if 'transcript_consequences' in info: 
            for i in info['transcript_consequences']:
                if most in  i['consequence_terms'] :
                    csqAllele=i['variant_allele']
                    freq_dict['csqAllele']=csqAllele
                    freq_dict['gene_id']=i['gene_id']
                    freq_dict['gene_symbol']=i['gene_symbol']
                else:
                    if 'regulatory_feature_consequences' in info: 
                        for r  in info['regulatory_feature_consequences']:
                            if most in  r['consequence_terms']: 
                                csqAllele=r['variant_allele']
                                freq_dict['csqAllele']=csqAllel
    return freq_dict
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def csqAlleleFeatures (csqAllele, altAllele, altAlleleCount, ploidy): 
	#csqAllele= cons allele  from vep  
	#altAlleleCount= integer, alternate allele count 
	# ploidy = number of allels at one locus 
	
	if csqAllele == altAllele: mycsqAlleleCount = altAlleleCount      ## Consequence allele is the alternate allele 
	else: mycsqAlleleCount = ploidy-altAlleleCount    ## Consequence allele is the reference allele;
	
	return [ csqAllele, mycsqAlleleCount]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def csqAlleleFeaturesLike (csqAllele, altAllele, altAlleleCount, GL ): 
	#csqAllele= cons allele  from vep  
	#altAlleleCount= integer, alternate allele count 
	#GL, string, comma separated vector of genotpe likelihoods
	#"""

	if csqAllele == altAllele: 
		mycsqAlleleCount = altAlleleCount      ## csqAllele counts  
	else: 
		mycsqAlleleCount = 2-altAlleleCount 
	Likl=GL.split(",")[altAlleleCount] ## csqAllele Likelihood 
	return [ csqAllele, mycsqAlleleCount, Likl]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def csqAlleleFeaturesMulti(genotype, csqAllele, refAllele, altAlleles, altAlleleCount, GL): 
	allAlleles=[refAllele]+ altAlleles.split(',')
	templist=altAlleleCount.split(',')
	allCounts=[2-sum( map(int, templist) )]+templist
	dCount=dict(zip(allAlleles, allCounts)) 
	if csqAllele in dCount:  
		csqCount=int(dCount[csqAllele]) 

		""" there is no itertool or numpy therfore I HAVE TO write the following...
		it MUST be rewritten to take the upperdiagonal 
		use numpy
		"""

		tempGenolist=[str(a)+str(b) for a in range (len(allAlleles)) for b in range (len(allAlleles)) ]
		"""This is equivalent to: 

		tempGenolist=[]
		for a in range (len(allAlleles)):		
			for b in range (len(allAlleles)): 
				tempGenolist.append(str(a) +str(b)) 
		"""
		if len(allAlleles)==2: genotypesIndexToRetain=[0,1,3]
		elif len(allAlleles)==3: genotypesIndexToRetain=[0,1,2,4,5,8]
		
		tempGenolistNoEquivalents=[i for i in tempGenolist if tempGenolist.index(i) in genotypesIndexToRetain] 
		dGL=dict(zip(tempGenolistNoEquivalents, GL.split(',')  ))	
		genotypeFormat=str(genotype[0])+str(genotype[2])
		likl=dGL[genotypeFormat]

		myout= [csqAllele, csqCount, likl]
	else: myout =False   ## no match between csq allele and alleles of teh vcf 	
	return myout

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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def CountCSQ_REF_ALT (csqAllele, refAllele, altAlleles, GTfields):
	"""csqAllele = type: string, consequence allele  from vep  
	refAllele = type: string, reference allele from variant calling
	altAlleles = type: string, comma-separated list of alternate allele from variant calling
	nbAploidSamples = type: int, number of total Alleles
	GTfields = type: list, linesplit[9:] (take only from 9° column to the end) ["0/0:.:.:.:.:.:.:.","1/1:16:0,16:0:0:16:628:-34.365,-4.81648,0"]
	hetGenotypes = type: int, heterozygosity in samples
	countPassLine = type: int, line the are PASS in frequency analisys	"""
	myres=False
	splitAltAlleles=altAlleles.split(",")
	allAlleles=[refAllele]+ splitAltAlleles
	mygstring=""; hetGenotypes=0; countPassLine=0;
	GTsplit=[i.split(":")[0] for i in GTfields]
	
	for idx, item in enumerate(GTsplit):
		if item !='./.':#if i is not "./.":
			mygstring+=item; 
			if item[0]!=item[2]: hetGenotypes+=1
	CountAlleles=[]
	for i in range(len(allAlleles)):
		#if str(i) in mygstring:
		CountAlleles.append(mygstring.count(str(i)))
	dAllele=dict(zip(allAlleles,CountAlleles))
		
	freqREF="{0:4.2f}".format(dAllele[refAllele])
        
        #### ALT freq 
	sumALT=0
	for i in splitAltAlleles:
		sumALT+=dAllele[i]
	freqALT="{0:4.2f}".format(sumALT)
	#### MAF freq
	MAF=min(freqALT,freqREF)
	#### CSQ freq
	if csqAllele in dAllele: 
		csqAllCount=dAllele[csqAllele]	
		freqCsq="{0:4.2f}".format(csqAllCount) 
			
	else: freqCsq='NA'

	#### define myres 
	myres= [freqCsq, freqREF, freqALT, MAF]
	return myres
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def AnnotateFreqCSQ_REF_ALT (csqAllele, refAllele, altAlleles, GTfields):
	"""csqAllele = type: string, consequence allele  from vep  
	refAllele = type: string, reference allele from variant calling
	altAlleles = type: string, comma-separated list of alternate allele from variant calling
	nbAploidSamples = type: int, number of total Alleles
	GTfields = type: list, linesplit[9:] (take only from 9° column to the end) ["0/0:.:.:.:.:.:.:.","1/1:16:0,16:0:0:16:628:-34.365,-4.81648,0"]
	hetGenotypes = type: int, heterozygosity in samples
	countPassLine = type: int, line the are PASS in frequency analisys	"""
	myres=False
	splitAltAlleles=altAlleles.split(",")
	allAlleles=[refAllele]+ splitAltAlleles
	mygstring=""; nbHaploidSamples=0; hetGenotypes=0; countPassLine=0;
	GTsplit=[i.split(":")[0] for i in GTfields]
	
	for idx, item in enumerate(GTsplit):
		if item !='./.':#if i is not "./.":
			mygstring+=item; nbHaploidSamples+=2
			if item[0]!=item[2]: hetGenotypes+=1
	CountAlleles=[]
	for i in range(len(allAlleles)):
		#if str(i) in mygstring:
		CountAlleles.append(mygstring.count(str(i)))
	dAllele=dict(zip(allAlleles,CountAlleles))
		
	if nbHaploidSamples!=0:
		#countPassLine+=1
		#### REF freq
		freqREF="{0:4.2f}".format(dAllele[refAllele]/nbHaploidSamples)
        
		#### ALT freq 
		sumALT=0
		for i in splitAltAlleles:
			sumALT+=dAllele[i]
		freqALT="{0:4.2f}".format(sumALT/nbHaploidSamples)
		
		#### MAF freq
		MAF=min(freqALT,freqREF)
		
		#### CSQ freq
		if csqAllele in dAllele: 
			csqAllCount=dAllele[csqAllele]	
			freqCsq="{0:4.2f}".format(csqAllCount/nbHaploidSamples) 
			# for obtain a '%' instead of integer
			# "{0:4.2f}%" : 
			# 4 = four number include the float ,
			# 2 = two float numbers,
			# f = for floating numbers; 
			# "{0:.0f}%" = if i want only '%' without float.
		
		else: freqCsq='NA'

		#### heterozygosity
		het= (hetGenotypes*2)/float(nbHaploidSamples)

		#### define myres 
		myres= [freqCsq, freqREF, freqALT, MAF, het]
	return myres
	#	else:
	#		myres= ['NA', freqREF, freqALT, MAF, hetGenotypes, countPassLine]
	#		return myres
	#else: 
	#	myres= ['NA', 'NA', 'NA', 'NA', 'NA', 'NA']
	#	return myres
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

