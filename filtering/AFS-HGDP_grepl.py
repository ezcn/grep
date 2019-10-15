import re 
import sys
sys.path.append('/mpba0/vcolonna/gianluca/TESI/MergedFreqScrip/greplib.py')
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
	args = parser.parse_args()
	#output = open(args.o,'w')
	#print(args) 


#############################################################


	##  READ VEP consequences rank ########
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
	random.seed(899)
	
	sampleToConsider=random.sample(EUR, 6)


##########~~~~~~~~~~~~~~  Loop of vcf lines 
	filemyres=open(args.o, 'w')  ## write a file.gz
	listOfErrors=[]
	dInfo={}
	column2retain=[]
	header=["Chr", "Pos","VariantClass", "CSQallele","CSQrank","Consequence","CSQfreq","REFfreq","ALTfreq","MAF","SampleID",'\n']
	filemyres.write("\t".join(map(str, header)))
	for line in gzip.open(args.f, 'r'):
		decodedLine=line.decode()  ## line.decode() is necessary to read encoded data using gzip in python3
		if re.match('##', decodedLine):
			if re.search("ID=CSQ" , decodedLine ):
				csqHeader=decodedLine.rstrip().split(":")[1].lstrip().rstrip("\">").split("|")		
			#filemyres.write(decodedLine)
		elif re.match('#', decodedLine):
			for ind in sampleToConsider:
				column2retain.append(decodedLine.split().index(ind))
		else:
			#print("this is a new line ") ## line split by  tab
			linesplit=decodedLine.rstrip().split()
			#print(linesplit)
			mychr=linesplit[0]; mypos=linesplit[1]; myref=linesplit[3]; myalt=linesplit[4] ## basic info

			##~~ split INFO field
			tempinfo=linesplit[7] 
			for i in tempinfo.split(";"):
				if re.search('=', i):
					temp=i.split("=")
					dInfo[temp[0]]=temp[1]
				#if re.search('NEGATIVE_TRAIN_SITE', i): pass
				#elif re.search('POSITIVE_TRAIN_SITE', i): pass

				else: pass 
					#dInfo[temp[0]]=temp[1]
					#print(dInfo)


			##~~ work on dInfo[CSQ]

			##~~ split for multiple consequences separated by ","
			multipleCsq=dInfo["CSQ"].split(",") 

			##~~ single consequence
			#print ('~~~  this is a consequence in a line ')
			CSQcount=0
			for mcsq in multipleCsq:
				CSQcount+=1
				myres=[]
				#myres+=[mychr, mypos]
				dCsq=dict(zip(csqHeader, mcsq.split("|") ))  #############    ALL VEP INFO 
				#print (dCsq) 
				#myres.append(dCsq['Existing_variation']) 
				
				#~~~~~~~~~~~  identify the allele with consequences
				#linesplit[7]=tempinfo # reset info field
				mycsqAllele=dCsq["Allele"]
				GTfields=[]
				GTfields+=[linesplit[column2retain[0]],linesplit[column2retain[1]],linesplit[column2retain[2]],linesplit[column2retain[3]],linesplit[column2retain[4]],linesplit[column2retain[5]]] 
				nbAploidSamples=len(GTfields)*2
				freqCSQ_REF_ALT=gp.AnnotateFreqCSQ_REF_ALT(mycsqAllele,myref, myalt, nbAploidSamples, GTfields)
				#print(freqCSQ_REF_ALT)
				#print (dCsq['Consequence'].split("&"))
				for c in dCsq['Consequence'].split("&"):
					#~~~~~~~~~~~~ assign severity score at the  most severe csq
					myindexes=[]
					myindexes.append(lSOTerm.index(c ))
					mostSevereCsq=lSOTerm[min(myindexes)]
					linesplit[7]=tempinfo	# reset info field 
					myres+=(linesplit[0],linesplit[1],dCsq["VARIANT_CLASS"],dCsq["Allele"],dSOTermRank[mostSevereCsq],c,freqCSQ_REF_ALT[0],freqCSQ_REF_ALT[1],freqCSQ_REF_ALT[2],freqCSQ_REF_ALT[3],sampleToConsider)
					filemyres.write("\t".join( map(str, myres ) ) )
					filemyres.write('\n')
					#linesplit[7]+=";CSQallele="
						#linesplit[7]+=dCsq["Allele"]
						#linesplit[7]+=";Consequence="
						#linesplit[7]+=c
						#linesplit[7]+=";CSQrank="
						#linesplit[7]+=dSOTermRank[mostSevereCsq]
						#linesplit[7]+=";IMPACT="
						#linesplit[7]+=dCsq["IMPACT"]
						#linesplit[7]+=";ExistingVariation="
						#linesplit[7]+=dCsq["Existing_variation"]
						#linesplit[7]+=";VariantClass="
						#linesplit[7]+=dCsq["VARIANT_CLASS"]
						#linesplit[7]+=";nbCSQ=" # for specify the number of CSQ allele
						#linesplit[7]+=str(CSQcount) # for specify the number of CSQ allele 
						#linesplit[7]+=";CSQfreq=%s" %(freqCSQ_REF_ALT[0])
						#linesplit[7]+=";REFfreq=%s" %(freqCSQ_REF_ALT[1])
						#linesplit[7]+=";ALTfreq=%s" %(freqCSQ_REF_ALT[2])
						#linesplit[7]+=";MAF=%s" %(freqCSQ_REF_ALT[3])
						#filemyres.write("\t".join(linesplit))
						#filemyres.write("\n")
	#fileToWrite=open(args.e, 'w')
	#for i in listOfErrors: fileToWrite.write( i )
	
if __name__ == "__main__":
	main()
