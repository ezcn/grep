import re 
import sys
#sys.path.append('/mpba0/vcolonna/gianluca/pythonScript/greplib.py')
sys.path.append('/home/enza/ezcngit/grep/filtering/greplib.py')
import greplib as gp
import argparse
import gzip
import random
import pandas as pd
#from __future__ import division

########################################################
def averagesFromFile(VEPannotatedVCF, column2retain , lSOTerm):
	"""
	VEPannotatedVCF= vcf annotated using VEP
	column2retain= list, index of samples to consider
		lSOTerm= list of Sequence Ontolgy terms that define variant consequences
	"""
	dInfo={};  dSOT={}; dImpact={}
	for line in gzip.open(args.f, 'r'):
		decodedLine=line.decode()  ## line.decode() is necessary to read encoded data using gzip in python3
		if re.match('##', decodedLine):
			if re.search('ID=CSQ' , decodedLine ): csqHeader=decodedLine.rstrip().split(':')[1].lstrip().rstrip('\'>').split('|')
		elif re.match('#', decodedLine):
			for ind in sampleToConsider:
				column2retain.append(decodedLine.split().index(ind))
		else:
			linesplit=decodedLine.rstrip().split()
			mychr=linesplit[0]; mypos=linesplit[1]; myref=linesplit[3]; myalt=linesplit[4] ## basic info
			##~~ split INFO field
			tempinfo=linesplit[7] 
			for i in tempinfo.split(';'):
				if re.search('=', i): # check if INFO fields has a value corresponding
					temp=i.split('=')
					dInfo[temp[0]]=temp[1]
					
				else: pass
					
			##~~ work on dInfo[CSQ]
			
			##~~ split for multiple consequences separated by ','
			multipleCsq=dInfo['CSQ'].split(',') 
			##~~ single consequence
			#print ('~~~  this is a consequence in a line ')
			for mcsq in multipleCsq:
				dCsq=dict(zip(csqHeader, mcsq.split('|') ))  #############    ALL VEP INFO 
					 
				#~~~~~~~~~~~  identify the allele with consequences
							
				mycsqAllele=dCsq['Allele']
				GTfields=[]
				for col in range(args.n): GTfields+=[linesplit[column2retain[col]]]
				#nbAploidSamples=len(GTfields)*2
				
				freqCSQ_REF_ALT=gp.AnnotateFreqCSQ_REF_ALT(mycsqAllele,myref, myalt, GTfields) # calculate allelic frequencies 
				
				for cons in dCsq['Consequence'].split('&'):
					#~~~~~~~~~~~~ assign severity score at the  most severe csq

                                                
					if freqCSQ_REF_ALT[0]!='NA':
						if not dCsq['Consequence'] in dSOT: dSOT[dCsq['Consequence']]=[0,0] #inizialize di dictionary with [counter, allele freq] if the key is not present 
						dSOT[dCsq['Consequence']][0]+=1 #add +1 to the counter  
						dSOT[dCsq['Consequence']][1]+=float(freqCSQ_REF_ALT[0]) #add the value of the consequence allele 

					else: listOfErrors.append((mychr, mypos,myref, myalt, dCsq["Allele"]) ) #to be printed in the error file to compare allele matching  
	
	CsqMeans=[mychr,freqCSQ_REF_ALT[4]]
	for vcsq in lSOTerm:
		if vcsq in dSOT: CsqMeans.append(dSOT[vel][1]/float(dSOT[vel][0]))
		else: CsqMeans.append('na')
	return CsqMeans 

########################################################
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', help='path to  input  file ',type=str,required=True)
	parser.add_argument('-v', help='path to table of vep consequences  ',type=str, required= True)	
	parser.add_argument('-o', help='path to output file  ',type=str, required= True)
	parser.add_argument('-e', help='path to error file',type=str,required=True)
	parser.add_argument('-m', help='path to metadata file',type=str,required=True)
	parser.add_argument('-c1', help='number of random cycle for pop1',type=int,required=False, default=1)
	parser.add_argument('-c2', help='number of random cycle for pop2',type=int,required=False, default=1)
	parser.add_argument('-s', help='seeds number',type=int,required=True)
	parser.add_argument("-n", help="number of cycles/replicates",type=int,required=True)
	parser.add_argument("-p1", help="population to analyze",type=str,required=True) #EUROPE
	parser.add_argument("-p2", help="population to analyze",type=str,required=True) #GREP
	args = parser.parse_args()
	#output = open(args.o,'w')
	#print(args) 


#############################################################

##### 0a. retrieve VEP ranking info   
	lSOTerm=gp.VepSOTermInfo(args.v)
			

##########~~~~~~~~~~~~~~ Read metadata

	Region=[]
	Sample=[]
	
	for line in open(args.m, 'r'):
		if re.match('sample', line):
			header= line.rstrip().split()
		else:
			other=line.rstrip().split()
			dMeta= dict(zip(header, other))
			Region.append(dMeta['region'])
			Sample.append(dMeta['sample'])

	dSampleRegion=dict(zip(Sample, Region))
	#### 
	pop1 = [key  for (key, value) in dSampleRegion.items() if value == args.p1]
	pop1sorted = sorted(EUR) ## needed for seed
	random.seed(args.s) ## need a sorted list of EUR
	pop2=[key  for (key, value) in dSampleRegion.items() if value == args.p2]
	pop2sorted = sorted(GREPtoretain)
	#sampleToConsider=random.sample(EUR, 6)

##########~~~~~~~~~~~~~~  Loop of vcf lines 

	#sys.stdout=open(args.o, 'w') 
	listOfErrors=[]
	#preHeader=['replicate','pop', 'chr','heterozigosity']
	#Header=preHeader+lSOTerm
	#print ("\t".join([i for i in Header]))

	##### 1. parse vcf to produce dVcf[mykey]=[myref, myqual, dFormat["GT"]]; mykey is  1:333333:/T (T is the alternate allele) "
	dVcf={}
	for line in gzip.open(args.f, 'r'):
		decodedLine=line.decode()  ## line.decode() is necessary to read encoded data using gzip in python3
		if not re.match('#', decodedLine):
			linesplit=decodedLine.rstrip().split()
			mychr=linesplit[0]; mypos=linesplit[1]; myref=linesplit[3]; myalt=linesplit[4]; myqual=float(linesplit[5]); altAlleles=myalt.split(",")
			tempformattitle=linesplit[8].split(":")
			tempformatcontent=linesplit[9].split(":")
			for altAl in altAlleles:
				mykey=mychr.lstrip("chr") + ":" + mypos + ":/" + altAl
				dVcf[mykey]=[myref, myalt, myqual, [gg.split(":")[0] for gg in linesplit[9:]] ] 

	##### 2. get VEP info for dVcf.keys()
	listOfErrors=[]
	#dVep={}
	for locusID in dVcf.keys(): 
		#print(locusID)
		dVepSearch=gp.getMostSevereCsqFromVep(locusID)
		#print("ho finito dVep")
		#prnt(dVepValue) 
		if type(dVepSearch) is not str: 
			dVcf[locusID]=dVepSearch
			if "csqAllele" in dVcf[locusID]:
				dVcf[locusID]["csqCount"]= gp.CountCSQ_REF_ALT(dVcf[locusID]["csqAllele"], dVcf[locusID][0], dVcf[locusID][1], [dVcf[locusID][3]]) [0]
			else:
				dVcf[locusID]["csqCount"] = np.nan
		else: 
			listOfErrors.append(locusID)
			
	df = pd.DataFrame(dVep).T
	print(dVcf) 
	print(df) 
"""	cycle=0
	while cycle < args.c1:
		cycle+=1
		#dInfo={};  dSOT={}; dImpact={}
		column2retain=[]
		sampleToConsider=random.sample(pop1sorted, args.n)
		pop=args.p1
		myresPop1=[cycle]
		myresPop1+=pop
		#vectorOfMeansPop1=averagesFromFile(args.f, column2retain ,  lSOTerm, )
		#myresPop1+=vectorOfMeansPop1
		#print ("\t".join(map(str, myresPop1) ))
		
		
	pop=args.p2
	myresPop2=[0]
	myresPop2+=pop
	vectorOfMeansPop2=averagesFromFile(args.f, pop2sorted ,  lSOTerm)
	myresPop2+=vectorOfMeansPop2
	print ("\t".join(map(str, myresPop2) ))
"""						
if __name__ == '__main__':
	main()
