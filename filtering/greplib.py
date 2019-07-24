#!/usr/bin/python3

def csqAlleleFeatures (csqAllele, altAllele, altAlleleCount, GL ): 
	#csqAllele= cons allele  from vep  
	#altAlleleCount= integer, alternate allele count 
	#GL, string, comma separated vector of genotpe likelihoods
	#"""

	if csqAllele == altAllele: mycsqAlleleCount = altAlleleCount      ## csqAllele counts  
	else: mycsqAlleleCount = 2-altAlleleCount 
	Likl=GL.split(",")[altAlleleCount] ## csqAllele Likelihood 
	return [ csqAllele, mycsqAlleleCount, Likl]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def checkFreq (listFreq, threshold): 
	rare=True
	if len(listFreq) >0: 
		for freq in listFreq:
			if freq > threshold:
				rare=False
				break
	else: rare="na"
	return(rare) 
