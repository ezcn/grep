import re 
import sys
sys.path.append('/mpba0/vcolonna/silvia/prioritiz/greplib.py')
import greplib as gp
import argparse
import gzip


########################################################
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", help="path to  input  file ",type=str,required=True)
	parser.add_argument("-i", help="threshold for SOTerm Impact  ", type=int,required=True)
	parser.add_argument("-r", help="threshold for rare variant definition ", type=float,required=True) 
	parser.add_argument("-v", help="path to table of vep consequences  ",type=str, required= True)   
	parser.add_argument("-o", help="path to output file  ",type=str, required= True)
	parser.add_argument("-e", help="path to error file",type=str,required=True)
	args = parser.parse_args()
	#output = open(args.o,'w')
	#print(args) 


	##  READ VEP consequences rank ########
	dRank={"HIGH":4, "LOW": 2, "MODERATE":3, "MODIFIER":1}
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

	#print (lSOTerm)
	lScores=list(reversed(range(len(lSOTerm)))) 
	#print (lScores) 
	dSOTermFineRank=dict(zip(lSOTerm, map(int, lScores) ))
	#print (dSOTermFineRank)


##########~~~~~~~~~~~~~~  Loop of vcf lines 
	filemyres=open(args.o, 'w')
	listOfErrors=[]
	dInfo={}
	header=["chr", "pos", "csqAllel", "csqAlleleCount", "GTLiklihood" , "ENSTID", "ImpactScore", "FineImpactScore", "rare","Embryo","GnomAD","CellCycle","DDD",'\n']
	filemyres.write("\t".join(map(str, header)))   

	#gzip.open(args.f, 'r')	
	for line in open(args.f, 'r'):
		if not re.match('#', line): 
			#print("this is a new line ") ## line split by  tab 
			linesplit=line.rstrip().split()
			
			mychr=linesplit[0]; mypos=linesplit[1]; myref=linesplit[3]; myalt=linesplit[4] ## basic info  
				
			if len(myalt.split(","))> 2:   #~~ excludes cases with more than two alt allele 
				pass 
			else:
	
				##~~ split INFO field
				tempinfo=linesplit[7] 
				for i in tempinfo.split(";"):  
					temp=i.split("=") 
					dInfo[temp[0]]=temp[1]

					
				##~~~ split FORMAT field
				tempformattitle=linesplit[8].split(":")
				tempformatcontent=linesplit[9].split(":")
				dFormat=dict(zip(tempformattitle, tempformatcontent))

				##~~ work on dInfo[CSQ]
				##~~ split for multiple consequences separated by ","
				multipleCsq=dInfo["CSQ"].split(",") 

				###~~ single consequence
				#print ('~~~  this is a consequence in a line ')
				for mcsq in multipleCsq:    
					myres=[]
					myres+=[mychr, mypos]
					dCsq=dict(zip(csqHeader, mcsq.split("|") ))  #############    ALL VEP INFO 
				
					#~~~~~~~~~~~  identify the allele with consequences
					mycsqAllele=dCsq["Allele"] 
					#~~~~~~~~~~~  csq allele features 
					featMultiOut=gp.csqAlleleFeaturesMulti( dFormat["GT"], mycsqAllele, myref, myalt, dInfo["AC"], dFormat["GL"] ) ## features of csqAll				
					if not featMultiOut:
						listOfErrors.append( '\t'.join([mychr, mypos, 'csq allele not matching', '\n']))  
						break 
					#print (featMultiOut) 
					else: myres+=featMultiOut
					
					#~~~~~~~~~~~ append ENSTID to myres	
					if  dCsq['Feature'] is not '':	myres.append(dCsq['Feature'])
					else: myres.append('na')
						
					#~~~~~~~~~~~~ assign severity score at the  most severe csq
					myindexes=[]
					for tl in dCsq['Consequence'].split("&"): 
						myindexes.append(lSOTerm.index(tl ))	
					mostSevereCsq=lSOTerm[min(myindexes)]
					myres.append( dSOTermRank[mostSevereCsq ]) ## score based on the impact of the consequence       	
					myres.append( dSOTermFineRank[mostSevereCsq ])
					
					#~~~~~~~~~~~~~~ find out if a variant is rare 
					freqlist = [float(x) for x in  [dCsq["AFR_AF"],dCsq["AMR_AF"],dCsq["EAS_AF"],dCsq["EUR_AF"],dCsq["SAS_AF"]] if x ] 
					rare = gp.checkFreq (freqlist, args.r) # check if it is a rare variant (af<thresh)  
					myres.append(rare)
									
					#~~~~~~~~~~ check if row have Embryo,CellCycle,DDD,GmomAD genes
					embryo=False ; DDD=False; cellcycle=False; gnomAD=False

					if re.search("annotation", line): embryo=True
					myres.append(embryo)
					if re.search("ANN3", line): gnomAD=True
					myres.append(gnomAD)
					if re.search("ANN2", line): cellcycle=True
					myres.append(cellcycle)
					if re.search("ANN1", line): DDD=True
					myres.append(DDD)
					

				#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	

				
					if dSOTermFineRank[mostSevereCsq ] > args.i and rare==True and cellcycle==True or dSOTermFineRank[mostSevereCsq ] > args.i and rare==True and embryo==True or dSOTermFineRank[mostSevereCsq ] > args.i and rare==True and DDD==True or dSOTermFineRank[mostSevereCsq ] > args.i and rare==True and gnomAD==True:	
						filemyres.write("\t".join( map(str, myres ) ) )
						filemyres.write('\n')

		else: 
			if re.search("ID=CSQ" ,line ): 
				csqHeader=line.rstrip().split(":")[1].lstrip().rstrip("\">").split("|")		
				#print (csqHeader)	

	fileToWrite=open(args.e, 'w')
	for i in listOfErrors: fileToWrite.write( i )
 


if __name__ == "__main__":
	main()
