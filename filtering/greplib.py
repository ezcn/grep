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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def csqAlleleFeaturesMulti(genotype, csqAllele, refAllele, altAlleles, altAlleleCount, GL): 
	allAlleles=[refAllele]+ altAlleles.split(',')
	templist=altAlleleCount.split(',')
	allCounts=[2-sum( map(int, templist) )]+templist
	dCount=dict(zip(allAlleles, allCounts)) 
	csqCount=dCount[csqAllele]

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

		
	return [csqAllele, csqCount, likl] 



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
