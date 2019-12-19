import requests, sys, json ,  argparse

""" 
USAGE: python3 vep.py -t <threshold below with variant is retained> -f <disgenet.test.out>

INPUT: an output file form disgenet.py or in general any file that contains rsid in the form 
      "rs1060501145",
(one or more rsid per line, rsdi must have quotes

OUTPUT: a single column text file with rsids of variants that are either not described in the genomic databses contained in VEP or for which there is at least one population described in VEP with allele frequency >= the threshold indicated in the command line   

"""
def compareFreq (dictionary, threshold): 
	variantIsCommon=0; retain=True 	
	for elem in dictionary:
		if len( [i  for i in dictionary[elem].values() if i>=threshold ]) > 0: variantIsCommon+=1 
	if variantIsCommon > 0 : retain =False  
	return retain 


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

def get_rsID (line): 
	linesp=line.rstrip().split("\"")
	myrs=[ i for i in linesp if "rs" in i] 
	return myrs[0]
 
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", help="path to  input  file ",type=str,required=True)
	parser.add_argument("-t", help="threshold for filtering allele frequencies ", type=float,required=True)
	args = parser.parse_args()	

	for line in open(args.f, 'r'):
		if "\"rs" in line: 
			rselem = get_rsID(line )
			mydict = getfreqfromVEP(rselem)
			if not any(mydict): print (rselem) #, mydict , "RARE")
			elif compareFreq(mydict, args.t ): print(rselem) #, mydict, "YAY, RARE"  ) 
			#else: print(mydict, "DISCARDED")

if __name__ == "__main__":
	main()
