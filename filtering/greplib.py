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
""" check if a variant is rare:  none of the populations in listfreq has frequency greater than threshold"""
	rare=True
	if len(listFreq) >0: 
		for freq in listFreq:
			if freq > threshold:
				rare=False
				break
	else: rare="NOB" #never observed 
	return(rare)

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

