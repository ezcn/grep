#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def Freq_CSQ_REF_ALT (csqAllele, refAllele, altAlleles, missing_data_format, genotypeslist):
	"""csqAllele = type: string, consequence allele  from vep  
	refAllele = type: string, reference allele from variant calling
	altAlleles = type: string, comma-separated list of alternate allele from variant calling
	nbAploidSamples : calculated 
	GTfields = type: list, list of genotypes ["0/0", "0/1", ".", "./."]
	hetGenotypes = type: int, heterozygosity in samples"""
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
		freqREF="{0:4.2f}".format(dAllele[refAllele]/nbHaploidSamples)
		freqAlt=[]
		for i in listAlt: 
			freqAltTemp="{0:4.2f}".format(dAllele[i]/nbHaploidSamples)
			freqAlt.append(freqAltTemp)
   		#~~~ CSQ 
		if csqAllele in dAllele:
			csqAllCount=dAllele[csqAllele]
			freqCsq="{0:4.2f}".format(csqAllCount/nbHaploidSamples)	
		else: freqCsq='NA'
		myres= [freqCsq, freqREF, freqAlt]
	return myres

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getInfoFromVepLocally (jsonWithVEPannotations)
	import json 
	vepInfo={}; vepInfoCommon={}
        #~~ filename is the json output of VEP runned locally 
        #~~ filename is parsed by line and by alternate allele  
	with open(jsonWithVEPannotations, 'r') as f:
                for line in f:
                        info=json.loads(line) #~~ make a dictionary out of a string 
                        #~~ derive the key chr:pos:/alt from the "input" element and process each alternate allele 
                        locusdata=info['input'].split(); altAlleles=locusdata[4].split(","); mychr=locusdata[0]; mypos=locusdata[1]
                        for altAl in altAlleles:
                                mykey=mychr.lstrip("chr") + ":" + mypos + ":/" + altAl
                                #print(mykey) 
                                most=info["most_severe_consequence"]
                                vepInfoCommon[mykey]={}
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
                                #~~ check if the most severe consequence is in an intergenic region 
                                elif 'intergenic_consequences' in info:
                                        for ic in info['intergenic_consequences']:
                                                if most in  ic['consequence_terms']:
                                                        vepInfo[(mykey, 'intergenic')]={}
                                                        icAllele=ic['variant_allele']
                                                        if icAllele ==altAl:
                                                                vepInfo[(mykey, 'intergenic') ]['csqAllele']=icAllele
                                                                vepInfo[(mykey, 'intergenic') ]['key']=mykey
                                                                vepInfo[(mykey, 'intergenic') ]['element_id']='intergenic'
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


    #print (json.dumps(info, indent=4 ) )  
    #print (json.dumps(vepInfo, indent=4) ) 
        return vepInfo, vepInfoCommon



