import re 
import sys
sys.path.append('../libraries')
import greplib as gp
import argparse
import gzip
import random
import pandas as pd
import numpy as np 
#from __future__ import division


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def replicatesResults (numberOfCycles, lsotermList, listOfIDs , numberOfIndividuals, vcfFilegz , dVepCommon) :    
	dRepPop={} # to store number,  mean, and  sd {consequence: ( {'n': 10, 'm': 11.8, 'sd': 2.4},  {'n': 20, 'm': 15.3, 'sd': 3.2})} 
	for  term in lsotermList: dRepPop[term]=[]
	cycle=0
	while cycle < numberOfCycles:
		cycle+=1
		tempRepPop={}
		for  term in lsotermList: tempRepPop[term]=[]
		#~~ select a random sample from pop1 at each cycle 
		column2retain=[]
		sampleToConsider=random.sample(listOfIDs,  numberOfIndividuals)
		#print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		#print( sampleToConsider) 	

		#~~ extract info from random sample from vcf 	
		for line in gzip.open(vcfFilegz,  'r'):
			decodedLine=line.decode()  ## line.decode() is necessary to read encoded data using gzip in python3
			if re.match ('#CHR', decodedLine): 
				for ind in sampleToConsider:
					column2retain.append(decodedLine.split().index(ind))
				#print( column2retain) 
			if not re.match('#', decodedLine):
				linesplit=decodedLine.rstrip().split()
				mychr=linesplit[0]; mypos=linesplit[1]; myref=linesplit[3]; myalt=linesplit[4]; altAlleles=myalt.split(",")
				#genotypesToConsider=[gg.split(":")[0] for gg in linesplit[9:] if linesplit[9:].index(gg)+9   in column2retain]
				genotypesToConsider=[]
				for indx in column2retain: genotypesToConsider.append(linesplit[indx].split(":")[0])

				#print('######################################')
				#print(linesplit[9:])
				#print (genotypesToConsider)
				for altAl in altAlleles:
					mykey=mychr.lstrip("chr") + ":" + mypos + ":/" + altAl
					most=dVepCommon[mykey]['most_severe_consequence']
					csqAllele=dVepCommon[mykey]['csqAllele']
					myfreq=gp.Freq_CSQ_REF_ALT (csqAllele, myref, myalt , "." ,genotypesToConsider)
					tempRepPop[most].append(float(myfreq[0])  ) 
					#print (tempRepPop) 
		for consType  in tempRepPop:
			templist=tempRepPop[consType] 
			n, m, sd = len(templist), np.mean(templist), np.std(templist) 
			dRepPop[consType].append({'n': n, 'm': m, 'sd': sd}) 

	return dRepPop					
	

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
	args = parser.parse_args()
	#output = open(args.o,'w')
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
	dRepPop1=replicatesResults (args.c1, lSOTerm, listofpops[0] , args.n, args.f , dVepCommon) 
	dRepPop2=replicatesResults (args.c2, lSOTerm, listofpops[1] , args.n, args.f , dVepCommon)

	print ('POP1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	print(dRepPop1)
	print ('POP2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	print(dRepPop2)

if __name__ == '__main__':
	main()
