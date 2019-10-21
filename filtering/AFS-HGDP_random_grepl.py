import re 
import sys
sys.path.append('/mpba0/vcolonna/gianluca/pythonScript/greplib.py')
import greplib as gp
import argparse
import gzip
import random
#from __future__ import division


########################################################
	
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", help="path to  input  file ",type=str,required=True)
	parser.add_argument("-v", help="path to table of vep consequences  ",type=str, required= True)	
	parser.add_argument("-o", help="path to output file  ",type=str, required= True)
	parser.add_argument("-e", help="path to error file",type=str,required=True)
	parser.add_argument("-m", help="path to metadata file",type=str,required=True)
	parser.add_argument("-c", help="number of random cycle",type=int,required=False, default=1)
	parser.add_argument("-s", help="seed's number",type=int,required=True)
	args = parser.parse_args()
	#output = open(args.o,'w')
	#print(args) 


#############################################################

##########~~~~~~~~~~~~~~ READ VEP consequences rank
	
	"""read external file with info on VEP consequences  """
	dRank={"HIGH":"HIGH", "LOW": "LOW", "MODERATE":"MODERATE", "MODIFIER":"MODIFIER"}
	dSOTermRank={}
	lSOTerm=[]  ### list of SOTerm ordered by severity

	countlinesCsq= True
	for csqLine in open(args.v, 'r'):
		if countlinesCsq:
			csqTitle=csqLine.rstrip().split('\t')
			countlinesCsq=False
		else:
			myRowList=csqLine.rstrip().split('\t')
			dCsq= dict(zip(csqTitle, myRowList ))
			dSOTermRank[dCsq['SO term']]=dRank[dCsq['IMPACT']]
			lSOTerm.append(myRowList[0])

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
	
	EUR = [key  for (key, value) in dSampleRegion.items() if value == 'EUROPE']

	EURsorted = sorted(EUR) ## needed for seed

	random.seed(args.s) ## need a sorted list of EUR

	#sampleToConsider=random.sample(EUR, 6)

##########~~~~~~~~~~~~~~  Loop of vcf lines 

	sys.stdout=open(args.o, 'w')  
	listOfErrors=[]
	print("Chr","\t", "Pos","\t","VariantClass","\t", "CSQallele","\t","CSQrank","\t","Consequence","\t","CSQfreq","\t","REFfreq","\t","ALTfreq","\t","MAF","\t","Cycle","\t","Population")
	cycle=0
	while cycle < args.c:
		cycle+=1
		myinput=gzip.open(args.f, 'r')
		dInfo={}
		column2retain=[]
		sampleToConsider=random.sample(EURsorted, 6)
		for line in myinput:
			decodedLine=line.decode()  ## line.decode() is necessary to read encoded data using gzip in python3
			if re.match('##', decodedLine):
				if re.search("ID=CSQ" , decodedLine ):
					csqHeader=decodedLine.rstrip().split(":")[1].lstrip().rstrip("\">").split("|")		
					
			elif re.match('#', decodedLine):
				for ind in sampleToConsider:
					column2retain.append(decodedLine.split().index(ind))
			else:
				#print("this is a new line ") ## line split by  tab
				linesplit=decodedLine.rstrip().split()
				mychr=linesplit[0]; mypos=linesplit[1]; myref=linesplit[3]; myalt=linesplit[4] ## basic info

				##~~ split INFO field
				tempinfo=linesplit[7] 
				for i in tempinfo.split(";"):
					if re.search('=', i): # check if INFO fields has a value corresponding
						temp=i.split("=")
						dInfo[temp[0]]=temp[1]
					
					else: pass 
					
				##~~ work on dInfo[CSQ]
				
				##~~ split for multiple consequences separated by ","
				multipleCsq=dInfo["CSQ"].split(",") 
				##~~ single consequence
				#print ('~~~  this is a consequence in a line ')
				CSQcount=0
				for mcsq in multipleCsq:
					CSQcount+=1
					dCsq=dict(zip(csqHeader, mcsq.split("|") ))  #############    ALL VEP INFO 
					 
					#~~~~~~~~~~~  identify the allele with consequences
							
					mycsqAllele=dCsq["Allele"]
					GTfields=[]
					for col in range(6): GTfields+=[linesplit[column2retain[col]]]
					nbAploidSamples=len(GTfields)*2
					freqCSQ_REF_ALT=gp.AnnotateFreqCSQ_REF_ALT(mycsqAllele,myref, myalt, nbAploidSamples, GTfields)

					for cons in dCsq['Consequence'].split("&"):
						#~~~~~~~~~~~~ assign severity score at the  most severe csq
						myindexes=[]
						myindexes.append(lSOTerm.index(cons))
						mostSevereCsq=lSOTerm[min(myindexes)]
						print(linesplit[0],"\t",linesplit[1],"\t",dCsq["VARIANT_CLASS"],"\t",dCsq["Allele"],"\t",dSOTermRank[mostSevereCsq],"\t",cons,"\t",freqCSQ_REF_ALT[0],"\t",freqCSQ_REF_ALT[1],"\t",freqCSQ_REF_ALT[2],"\t",freqCSQ_REF_ALT[3],"\t",cycle,"\t","HGDP")


	#fileToWrite=open(args.e, 'w')
	#for i in listOfErrors: fileToWrite.write( i )

	
if __name__ == "__main__":
	main()
