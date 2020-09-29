import sys, json 

jsonWithVEPannotations=sys.argv[1] 
#thresholdRareAllele=sys.argv[2]

def getInfoFromVepLocally (jsonWithVEPannotations):
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
				most=vepInfo[mykey]=info["most_severe_consequence"]
				#~~  check if the most sequence is in a transcript
				if 'transcript_consequences' in info:
					for tc in info['transcript_consequences']:
						if most in  tc['consequence_terms'] :
							tcTranscript=info['transcript_consequences']['transcript_id']
							vepInfo[(mykey, tcTranscript)]={}
							tcAllele=tc['variant_allele']
							if tcAllele ==altAl :
						    		vepInfo[(mykey, tcTranscript)]['csqAllele']=tcAllele
						    		vepInfo[(mykey, tcTranscript)]['gene_id']=tc['gene_id']
						    		vepInfo[(mykey, tcTranscript)]['gene_symbol']= tc['gene_symbol']
						    		vepInfo[(mykey, tcTranscript)]['impact']=tc['impact']
				#~~ check if the most severe consequence is in a regulatory feature  
				elif 'regulatory_feature_consequences' in info:
					for rf  in info['regulatory_feature_consequences']:
						if most in  rf['consequence_terms']:
							rfRegulatory=info['regulatory_feature_id']
							vepInfo[(mykey, rfRegulatory)]={}
							rfAllele=rf['variant_allele']
							if rfAllele ==altAl :
								vepInfo[(mykey, rfRegulatory)]['csqAllele']=rfAllele
								vepInfo[(mykey, rfRegulatory)]['impact']=rf['impact']
				#~~ check if the most severe consequence is in an intergenic region 
				elif 'intergenic_consequences' in info:
					for ic in info['intergenic_consequences']:
						if most in  ic['consequence_terms']:
							vepInfo[(mykey, 'intergenic')]={}
							icAllele=ic['variant_allele']
							if icAllele ==altAl: 	
								vepInfo[(mykey, 'intergenic') ]['csqAllele']=icAllele

				#~~ retrive info on rsID, starting position, frequencies for the csqAllele found in the previous code  
				if "colocated_variants" in info:
					infoCV = info["colocated_variants"][0]
					vepInfoCommon[mykey]={}
					#~~ adding info for starting position, this is needed for merging with CADD score.
					vepInfoCommon[mykey]["start"] = infoCV["start"]

					#~~ check if rsid is present 
					if "id" in infoCV: vepInfoCommon[mykey]["id"] = infoCV["id"]
					
					#~~ check if variant is rare according to threshold given as argument and 
					#if 'frequencies' in infoCV:
					#	frqInfoDic=infoCV["frequencies" ]
					#	for i in frqInfoDic:
					 #   	if 'csqAllele' in vepInfoCommon[mykey]:
					#		if vepInfo[mykey]['csqAllele'] in frqInfoDic:
					#	    	to_parse = list(frqInfoDic[vepInfo[mykey]['csqAllele']].values())
					#	    	#print(to_parse)
					#	    	vepInfo[mykey]['rare']= checkFreq (to_parse, thresholdRareAllele )
					#				    #else: vepInfo[mykey]['rare']='na'


    #print (json.dumps(info, indent=4 ) )  
    #print (json.dumps(vepInfo, indent=4) ) 
	return vepInfo, vepInfoCommon



cicci= getInfoFromVepLocally(jsonWithVEPannotations) 
print(cicci) 
