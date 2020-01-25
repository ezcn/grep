import json, pprint 
filename="/lustre/home/enza/vepJson/AS006/AS006.chr1.vep.json"
filename="/lustrehome/enza/ezcngit/grep/filtering/cicci"


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

vepInfo={}
#~~ filename is the json output of VEP runned locally 
#~~ filename is parsed by line and by alternate allele  
with open(filename, 'r') as f:  
	for line in f:
		info=json.loads(line) #~~ make a dictionary out of a string 
		#print(line_d['input']) 
		#~~ derive the key chr:pos:/alt from the "input" element and process each alternate allele 
		locusdata=info['input'].split(); altAlleles=locusdata[4].split(","); mychr=locusdata[0]; mypos=locusdata[1]
		for altAl in altAlleles:
			mykey=mychr.lstrip("chr") + ":" + mypos + ":/" + altAl
			vepInfo[mykey]={}	

			#~~ check if the most severe consequence is in the line else retrun an empty vepInfo[mykey] 
			if "most_severe_consequence" in info:
				vepInfo[mykey]["most_severe_consequence"]=info["most_severe_consequence"]
				most=vepInfo[mykey]["most_severe_consequence"]
				#~~  check if the most sequence is in a transcript else in a regulatory feature 
				if 'transcript_consequences' in info:
					ReportedGeneId=[] ; ReportedGeneSymbol=[]
					for tc in info['transcript_consequences']:
						if most in  tc['consequence_terms'] :
							csqAllele=tc['variant_allele']
							if csqAllele ==altAl : 
								print (most, csqAllele, tc["transcript_id"], tc['gene_id'])
								vepInfo[mykey]['csqAllele']=csqAllele
								if not tc['gene_id'] in ReportedGeneId: ReportedGeneId.append(tc['gene_id'])
								if not tc['gene_symbol'] in ReportedGeneSymbol:  ReportedGeneSymbol.append(tc['gene_symbol']) 
					vepInfo[mykey]['gene_id']=ReportedGeneId
					vepInfo[mykey]['gene_symbol']=ReportedGeneSymbol

				elif 'regulatory_feature_consequences' in info: 
					for rf  in info['regulatory_feature_consequences']:
						if most in  rf['consequence_terms']:
							csqAllele=rf['variant_allele']
							if csqAllele ==altAl :
								vepInfo[mykey]['csqAllele']=csqAllel

			#~~ retrive info on rsID, starting position, frequencies
			if "colocated_variants" in info:
				id_search = info["colocated_variants"][0]
				if "id" in id_search:
			    		vepInfo[mykey]["id"] = id_search["id"]
				adding info for starting position, this is needed for merging with CADD score.
				if "start" in id_search: 
					vepInfo[mykey]["starting_position"] = id_search["start"]   
				#~~ check if variant is rare according to threshold given as argument and report True or False in vepInfo[mykey]
				if 'frequencies' in id_search:
					frqInfoDic=id_search["frequencies" ]
					for i in frqInfoDic: 
						if vepInfo[mykey]['csqAllele'] in frqInfoDic: 
							to_parse = list(frqInfoDic[vepInfo[mykey]['csqAllele']].values())
							#print(to_parse)
							vepInfo[mykey]['rare']= checkFreq (to_parse, 0.01) 		
							
 	

print (json.dumps(info, indent=4 ) )  
print (json.dumps(vepInfo, indent=4) ) 
