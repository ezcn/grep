import re 
import sys
sys.path.append('/mpba0/vcolonna/gianluca/TESI/MergedFreqScrip/greplib.py')
import greplib as gp
import argparse
import gzip
#from __future__ import division


########################################################
	
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", help="path to  input  file ",type=str,required=True)
	parser.add_argument("-o", help="path to output file  ",type=str, required= True)
	parser.add_argument("-e", help="path to error file",type=str,required=True)
	args = parser.parse_args()
	#output = open(args.o,'w')
	#print(args) 


##########~~~~~~~~~~~~~~  Loop of vcf lines 
	filemyres=open(args.o, 'w')
	listOfErrors=[]
	dInfo={}

	for line in gzip.open(args.f, 'r'):
		decodedLine=line.decode()  ## line.decode() is necessary to read encoded data using gzip in python3
		if re.match('##', decodedLine):
			if re.search("ID=CSQ" , decodedLine ):
				csqHeader=decodedLine.rstrip().split(":")[1].lstrip().rstrip("\">").split("|")		
			filemyres.write(decodedLine)
			
		elif re.match('#', decodedLine): 
			filemyres.write('##INFO=<ID=CSQfreq,Number=A,Type=Float,Description="Frequency of CSQ allele in set of samples">')
			filemyres.write("\n")
			filemyres.write('##INFO=<ID=REFfreq,Number=A,Type=Float,Description="Frequency of REF allele in set of samples">')
			filemyres.write("\n")
			filemyres.write('##INFO=<ID=ALTfreq,Number=A,Type=Float,Description="Frequency of ALT allele in set of samples">')
			filemyres.write("\n")
			filemyres.write('##INFO=<ID=nbCSQ,Number=A,Type=Integer,Description="1st, 2nd .... CSQ allele">')
			filemyres.write("\n")
			filemyres.write(decodedLine)
		else:
			#print("this is a new line ") ## line split by  tab
			linesplit=decodedLine.rstrip().split()
			#print(linesplit)
			mychr=linesplit[0]; mypos=linesplit[1]; myref=linesplit[3]; myalt=linesplit[4] ## basic info  

			#nbOfAltAlleles=len(myalt.split(","))			

			#if nbOfAltAlleles> 2:   #~~ excludes cases with more than two alt allele, want to add the excluded in the error output 
				#listOfErrors+=( '\t'.join([mychr, mypos, 'more than two alternate alleles', nbOfAltAlleles ,  '\n']))
				#pass 
			#else:
			##~~ split INFO field
			tempinfo=linesplit[7] 
			for i in tempinfo.split(";"):  
				temp=i.split("=") 
				dInfo[temp[0]]=temp[1]


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
				linesplit[7]=tempinfo # reset info field
				mycsqAllele=dCsq["Allele"]
				GTfields=linesplit[9:] 
				nbAploidSamples=len(GTfields)*2
				freqCSQ_REF_ALT=gp.AnnotateFreqCSQ_REF_ALT(mycsqAllele,myref, myalt, nbAploidSamples, GTfields)
				#print(freqCSQ_REF_ALT)
				linesplit[7]+=";"
				linesplit[7]+="nbCSQ=" # for specify the number of CSQ allele
				linesplit[7]+=str(CSQcount) # for specify the number of CSQ allele 
				linesplit[7]+=";CSQfreq=%s" %(freqCSQ_REF_ALT[0])
				linesplit[7]+=";REFfreq=%s" %(freqCSQ_REF_ALT[1])
				linesplit[7]+=";ALTfreq=%s" %(freqCSQ_REF_ALT[2])
				#CSQfreq="CSQfreq=%s" %(freqCSQ_REF_ALT[0])
				#REFfreq="REFfreq=%s" %(freqCSQ_REF_ALT[1])
				#ALTfreq="ALTfreq=%s" %(freqCSQ_REF_ALT[2])
				#allfreq=CSQfreq,REFfreq,ALTfreq
				#freqField=";".join(map(str,allfreq))
				#rowWithFreq=linesplit[7],freqField
				#row2print=";".join(map(str,rowWithFreq))
				#linesplit[7]=row2print
				filemyres.write("\t".join(linesplit))
				filemyres.write("\n")
	#fileToWrite=open(args.e, 'w')
	#for i in listOfErrors: fileToWrite.write( i )
	
if __name__ == "__main__":
	main()
