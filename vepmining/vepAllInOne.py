import requests, sys, json ,  argparse

""" 
send a request to https://rest.ensembl.org/documentation/info/vep_id_post to retrive information about allele frequencies in the general population (if available) for a list of rsid  

USAGE: python3 vepAllInOne.py -f <single column of rsid input file >   -t < allele frequency threshold below with variant is retained > 

INPUT: single column of rsid input file

OUTPUT: in the standard out, a single column of rsids of variants that are either not described in the genomic databses contained in VEP or for which there is at least one population described in VEP with allele frequency >= the threshold indicated in the command line   

"""

def compareFreq (dictionary, threshold): 
	variantIsCommon=0; retain=True 	
	for elem in dictionary:
		if len( [i  for i in dictionary[elem].values() if i>=threshold ]) > 0: variantIsCommon+=1 
	if variantIsCommon > 0 : retain =False  
	return retain 


def getfreqfromVEPbulck (listOfrsid):
	listOfDict=[]
	server = "https://rest.ensembl.org"
	ext = "/vep/human/id"
	headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
	strOfDataDictionary=json.dumps({"ids" : listOfrsid })
	res = requests.post(server+ext, headers=headers, data=strOfDataDictionary, params={"fields": "colocated_variants"})
	decoded = res.json()  # a python dictionary
	for elem in decoded:
		rsid=elem["input"]
		freq_dict={}
		for var in elem["colocated_variants"]:
			#rsid=var["id"]
			if 'frequencies' in var: freq_dict=var["frequencies"] 
		listOfDict.append((rsid, freq_dict))

	return listOfDict
	#return decoded

def getfreqfromVEP (rsid):
	freq_dict={} 
	server = "https://rest.ensembl.org"
	ext = "/vep/human/id/"+ rsid + "?"
	r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
	if not r.ok: 
		r.raise_for_status()
		sys.exit() 
	decoded = r.json()
	for var in decoded[0]["colocated_variants"] :
		if "frequencies" in var: freq_dict=var["frequencies"]
	return freq_dict 

 
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", help="path to  input  file ",type=str,required=True)
	parser.add_argument("-t", help="threshold for filtering allele frequencies ", type=float,required=True)
	args = parser.parse_args()		

	with open(args.f) as myf:
  		allrsId = myf.readlines()	
	#thislist=["rs12720452", "rs1060501130"]
	maxNumbRsIdForRequest=200 # there is a limit in the VEP rest  

	while len(allrsId) > 0 :
		sublist=allrsId[0:maxNumbRsIdForRequest]
		allrsId=allrsId[maxNumbRsIdForRequest:]
		result = getfreqfromVEPbulck(sublist) 
		for rsinfo in result: 
			rsdict=rsinfo[1]; rsname=rsinfo[0]
			if not any(rsdict): print (rsname) #, rsdict, "RARE")
			elif compareFreq(rsdict, args.t ): print(rsname) #, rsdict, "YAY, RARE"  )
			#else: print(rsname, rsdict, "DISCARDED")

		#print(json.dumps(result, indent=2))
		#print (result) 

if __name__ == "__main__":
	main()
