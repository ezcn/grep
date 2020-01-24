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
with open(filename, 'r') as f:
	for line in f:
		#print('#####')
		info=json.loads(line)
		#print(line_d['input']) 
		locusdata=info['input'].split(); altAlleles=locusdata[4].split(",")
		for altAl in altAlleles:
			mykey=locusdata[0].lstrip("chr") + ":" + locusdata[1] + ":/" + altAl
			vepInfo[mykey]={}	

			if "most_severe_consequence" in info:
				vepInfo[mykey]["most_severe_consequence"]=info["most_severe_consequence"]
				most=vepInfo[mykey]["most_severe_consequence"]
				if 'transcript_consequences' in info:
					ReportedGeneId=[] ; ReportedGeneSymbol=[]
					for i in info['transcript_consequences']:
						if most in  i['consequence_terms'] :
							csqAllele=i['variant_allele']
							if csqAllele ==altAl : 
								print (most, csqAllele, i["transcript_id"], i['gene_id'])
								vepInfo[mykey]['csqAllele']=csqAllele
								if not i['gene_id'] in ReportedGeneId: ReportedGeneId.append(i['gene_id'])
								if not i['gene_symbol'] in ReportedGeneSymbol:  ReportedGeneSymbol.append(i['gene_symbol']) 
					vepInfo[mykey]['gene_id']=ReportedGeneId
					vepInfo[mykey]['gene_symbol']=ReportedGeneSymbol

				elif 'regulatory_feature_consequences' in info: 
					for r  in info['regulatory_feature_consequences']:
						if most in  r['consequence_terms']:
							csqAllele=r['variant_allele']
							if csqAllele ==altAl :
								vepInfo[mykey]['csqAllele']=csqAllel
			if "colocated_variants" in info:
				id_search = info["colocated_variants"][0]
				if "id" in id_search:
			    		vepInfo[mykey]["id"] = id_search["id"]
				# adding info for starting position, this is needed for merging with CADD score.
				#if "start" in id_search: 
				#	vepInfo[mykey]["starting_position"] = id_search["start"]   
				if 'frequencies' in id_search:
					frqInfoDic=id_search["frequencies" ]
					for i in frqInfoDic: 
						if vepInfo[mykey]['csqAllele'] in frqInfoDic: 
							
							#to_parse = list(id_search["frequencies"].values())[0]
							to_parse = list(frqInfoDic[vepInfo[mykey]['csqAllele']].values())#[0]
							#print(to_parse)
							vepInfo[mykey]['rare']= checkFreq (to_parse, 0.01) 		
							
 	
					#for var in to_parse.keys():
					#	vepInfo[mykey][var]=to_parse[var]



			

#print (json.dumps(info, indent=4 ) )  
print (json.dumps(vepInfo, indent=4) ) 
#datastore = json.dumps(filecontent)
#for i in mydata: 
#	print (i ) 
