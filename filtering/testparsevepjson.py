import json, pprint 
filename="/lustre/home/enza/vepJson/AS006/AS006.chr1.vep.json"
filename="/lustrehome/enza/ezcngit/grep/filtering/cicci"

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
			if "colocated_variants" in info:
				id_search = info["colocated_variants"][0]
				if "id" in id_search:
			    		vepInfo[mykey]["id"] = id_search["id"]
				# adding info for starting position, this is needed for merging with CADD score.
				#if "start" in id_search: 
				#	vepInfo[mykey]["starting_position"] = id_search["start"]   
				if 'frequencies' in id_search:
					to_parse = list(id_search["frequencies"].values())[0]
					for var in to_parse.keys():
						vepInfo[mykey][var]=to_parse[var]

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

				if 'regulatory_feature_consequences' in info: 
					for r  in info['regulatory_feature_consequences']:
						if most in  r['consequence_terms']:
							csqAllele=r['variant_allele']
							if csqAllele ==altAl :
								vepInfo[mykey]['csqAllele']=csqAllel


#print (json.dumps(info, indent=4 ) )  
#mydata=json.loads(json.dumps(filecontent))
print (json.dumps(vepInfo, indent=4) ) 
#datastore = json.dumps(filecontent)
#for i in mydata: 
#	print (i ) 
