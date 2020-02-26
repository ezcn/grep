import re 
import sys
#sys.path.append('../libraries')
import grepGianluca as gp
import argparse
import gzip
import random
import pandas as pd
import numpy as np 
#from __future__ import division


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def replicatesResults (numberOfCycles, lsotermList, listOfIDs , numberOfIndividuals, vcfFilegz , dVepCommon, variantClass, rareTreshold) :    
	dRepPop={} # to store number,  mean, and  sd {consequence: ( {'n': 10, 'm': 11.8, 'sd': 2.4},  {'n': 20, 'm': 15.3, 'sd': 3.2})} 
	for  term in lsotermList: dRepPop[term]=[]
	cycle=0
	while cycle < numberOfCycles:
		cycle+=1
		tempRepPop={}
		for  term in lsotermList: tempRepPop[term]=[]
		#~~ select a random sample from listOfIDs at each cycle 
		column2retain=[]
		sampleToConsider=random.sample(listOfIDs,  numberOfIndividuals)

		#~~ extract info from random sample from vcf 	
		csqnotfound=0
		listnotfound=[]
		myf=gzip.open(vcfFilegz,  'r')
		for line in myf:
			decodedLine=line.decode()  ## line.decode() is necessary to read encoded data using gzip in python3
			if re.match ('#CHR', decodedLine):
				for ind in sampleToConsider:
					column2retain.append(decodedLine.split().index(ind))
			else: pass  
			
		for mykey in dVepCommon: 
			#~~~~~~~~ include the variant  based on variant class 
			if dVepCommon[mykey]['variant_class']==variantClass:

				#~~~~~~~  include the variant based on its frequency in gnomeAD and 1000G	
				if 'frequencies' in dVepCommon[mykey]: 
					for cAll in dVepCommon[mykey]['frequencies']: 
						listFreq=dVepCommon[mykey]['frequencies'][cAll].values()
						rare=gp.checkFreq (listFreq,  rareTreshold )
				#else:
				#	dVepCommon[celem]['rare']="NOB"
					if rare:  
						#~~~~~~~~  retain only genotypes of selected individuals 
						genotypesToConsider=[]
						for indx in column2retain: genotypesToConsider.append(dVepCommon[mykey]['genotypes'][indx-9])   #linesplit[indx].split(":")[0])
	
						most=dVepCommon[mykey]['most_severe_consequence']
						myfreq=gp.Freq_CSQ_REF_ALT ( mykey.split("/")[-1], dVepCommon[mykey]['ref'], mykey.split("/")[-1], "." ,genotypesToConsider)
			
						#~~~~~~~  
						if float(myfreq[0])>0: tempRepPop[most].append(float(myfreq[0]) ) 
					else: 
						listnotfound.append(dVepCommon[mykey]['frequencies'])

		for consType  in tempRepPop:
			templist=tempRepPop[consType] 
			n, m, sd = len(templist), np.mean(templist), np.std(templist) 
			dRepPop[consType].append({'n': n, 'm': m, 'sd': sd}) 

	return dRepPop, csqnotfound, listnotfound, tempRepPop					
	

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def MakeListOfIDs(populationMetadata, pop1name, pop2name): 
	""" at least two columns file with mandatory header "region" and "sample" """
	Region=[]; Sample=[]
	for line in open(populationMetadata, 'r'):
		if re.match('sample', line): header= line.rstrip().split()
		else:
			other=line.rstrip().split()
			dMeta= dict(zip(header, other))
			Region.append(dMeta['region'])
			Sample.append(dMeta['sample'])

	dSampleRegion=dict(zip(Sample, Region))
	pop1 = [key  for (key, value) in dSampleRegion.items() if value == pop1name]
	pop1sorted = sorted(pop1) ## needed for seed
	pop2=[key  for (key, value) in dSampleRegion.items() if value == pop2name]
	pop2sorted = sorted(pop2)
	return pop1sorted, pop2sorted 


#######################################################
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-j", help="path to  json  file ",type=str,required=True)
	parser.add_argument('-f', help='path to  vcf file (merged pop1 and pop2 ) ',type=str,required=True)
	parser.add_argument('-v', help='path to table of vep consequences  ',type=str, required= True)
	parser.add_argument('-o', help='path to output file  ',type=str, required= True)
	parser.add_argument('-e', help='path to error file',type=str,required=True)
	parser.add_argument('-m', help='path to metadata file',type=str,required=True)
	parser.add_argument('-c1', help='number of random cycle for pop1',type=int,required=False, default=1)
	parser.add_argument('-c2', help='number of random cycle for pop2',type=int,required=False, default=1)
	parser.add_argument('-s', help='seeds number',type=int,required=True)
	parser.add_argument("-n", help="number of individuals in pop1 and pop2",type=int,required=True)
	parser.add_argument("-p1", help="reference population fro standardization",type=str,required=True) #EUROPE
	parser.add_argument("-p2", help="test population ",type=str,required=True) #GREP
	parser.add_argument("-so", help="variant class as SO term. possible values:  SNV , insertion, deletion, indel, substitution , sequence_alteration ",type=str,required=True)
	parser.add_argument('-r', help='threshold for rare variant ',type=float,required=True)
	args = parser.parse_args()
	sys.stdout = open(args.o,'w')
	#print(args) 

	#############################################################

	#########~~~~~~~~~~~~ 0. retrieve VEP ranking info   
	lSOTerm=gp.VepSOTermInfo(args.v)

	#########~~~~~~~~~~~~ 1. get VEP info from local json file 
	dV = gp.getInfoFromVepLocally (args.j )   #yelds two dictionaries 
	dVepTrans=dV[0]; dVepCommon=dV[1]
		
	##########~~~~~~~~~~~ 2.  Read populations metadata
	listofpops=MakeListOfIDs(args.m, args.p1, args.p2)

	##########~~~~~~~~~~~ 3. process Vcfs 
	random.seed(args.s) # needed for ??? 
	dRepPop1=replicatesResults (args.c1, lSOTerm, listofpops[0] , args.n, args.f , dVepCommon, args.so, args.r ) 
	dRepPop2=replicatesResults (args.c2, lSOTerm, listofpops[1] , args.n, args.f , dVepCommon, args.so, args.r )

	#print ('POP1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	#print(dRepPop1[1])
	#print(dRepPop1[2]['upstream_gene_variant'])
	#print(dRepPop1[0])

	#print ('POP2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	#print(dRepPop2[1])
	#print(dRepPop2[2]['upstream_gene_variant'])
	#print(dRepPop2[0])
	#print (dRepPop2[3])
	xfile=open(args.e, "w")
	xfile.write(str(dRepPop1[2])  )

	musdHGDP={}
	for elem in dRepPop1[0]: 
		musdHGDP[elem]=gp.combineMeanSD( dRepPop1[0][elem])
	#print (musdHGDP)
	print('type' , 'n', 'muZFreq', 'sdZFreq', 'ci90ZFreq', 'ci95ZFreq', 'ci99ZFreq')
	for elem2 in dRepPop2[3]: 
		zfreq=[(x- musdHGDP[elem2][0])/musdHGDP[elem2][1] for x in  dRepPop2[3][elem2]]
		muZFreq=np.mean(zfreq)
		sdZFreq=np.std(zfreq)
		ci90ZFreq=1.645*(sdZFreq/np.sqrt(len(zfreq)))
		ci95ZFreq=1.960*(sdZFreq/np.sqrt(len(zfreq)))
		ci99ZFreq=2.576*(sdZFreq/np.sqrt(len(zfreq))) 
		#print('####')
		print (elem2 , len(zfreq), muZFreq, sdZFreq, ci90ZFreq, ci95ZFreq, ci99ZFreq) 
		#print ( zfreq)
	


if __name__ == '__main__':
	main()
