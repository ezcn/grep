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
	parser.add_argument("-c", help="consequence allele count ",type=int,required=False, default=0)
	parser.add_argument("-w", help="path to  weight  file ",type=str,required=True)
	args = parser.parse_args()
	#output = open(args.o,'w')
	#print(args) 
	

	## READ weights 
	dWeig={}
	for wline in open (args.w):
		w=wline.rstrip().split() 
		dWeig[w[1]]=int(w[2])
	#print (dWeights)
        ##  READ VEP consequences rank ########
	"""read external file with info on VEP consequences  """
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
	header=["chr", "pos", "Existing_variation", "quality", "csqAllel", "csqAlleleCount", "GTLiklihood" , "ENSTID", "ImpactScore", "FineImpactScore", "rare","Embryo","GnomAD","CellCycle","DDD","miscarriage","lethal", "essential", "gpscore", "wAvGpScore"]
	filemyres.write("\t".join(map(str, header)))
	filemyres.write("\n")   
	
	for line in gzip.open(args.f, 'r'):
		decodedLine=line.decode()  ## line.decode() is necessary to read encoded data using gzip in python3
		if re.match('#', decodedLine):
			if re.search("ID=CSQ" , decodedLine ):
				csqHeader=decodedLine.rstrip().split(":")[1].lstrip().rstrip("\">").split("|")		
			#print (csqHeader)	

		else:
			#print("this is a new line ") ## line split by  tab
			linesplit=decodedLine.rstrip().split()
			
			mychr=linesplit[0]; mypos=linesplit[1]; myref=linesplit[3]; myalt=linesplit[4]; myqual=linesplit[5] ## basic info  

			nbOfAltAlleles=len(myalt.split(","))			

			if nbOfAltAlleles> 2:   #~~ excludes cases with more than two alt allele, want to add the excluded in the error output 
				listOfErrors.append( '\t'.join([mychr, mypos, 'more than two alternate alleles', nbOfAltAlleles ,  '\n']))
				pass 
			else:
				##~~ split INFO field
				tempinfo=linesplit[7] 
				for i in tempinfo.split(";"):  
					temp=i.split("=") 
					dInfo[temp[0]]=temp[1]
				#print (dInfo) 	
				##~~~ split FORMAT field
				tempformattitle=linesplit[8].split(":")
				tempformatcontent=linesplit[9].split(":")
				dFormat=dict(zip(tempformattitle, tempformatcontent))

				##~~ work on dInfo[CSQ]
				##~~ split for multiple consequences separated by ","
				multipleCsq=dInfo["CSQ"].split(",") 

				##~~ single consequence
				#print ('~~~  this is a consequence in a line ')
				for mcsq in multipleCsq:    
					myres=[]; gpScore=0; gpweig=0
					myres+=[mychr, mypos]
					dCsq=dict(zip(csqHeader, mcsq.split("|") ))  #############    ALL VEP INFO 
					#print (dCsq) 
					myres.append(dCsq['Existing_variation']) # Existing_variation = identificativo 'rs'
					#~~~~~~~~~~~ append quality
					myres.append(myqual)										
					#~~~~~~~~~~~  identify the allele with consequences
					mycsqAllele=dCsq["Allele"] 
					#~~~~~~~~~~~  csq allele features : number of allele with consequences and genotype likelihood  
					featMultiOut=gp.csqAlleleFeaturesMulti( dFormat["GT"], mycsqAllele, myref, myalt, dInfo["AC"], dFormat["GL"] ) ## features of csqAll				
					if not featMultiOut:
						listOfErrors.append( '\t'.join([mychr, mypos, 'csq allele not matching', '\n']))  
						break 
					#print (featMultiOut) 
					else: 
						myres+=featMultiOut; 
						csqAllCount=int(featMultiOut[1])
						gpScore+=featMultiOut[1] * dWeig['wCAC']
						gpweig+=dWeig['wCAC']	
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
					gpScore+=dSOTermFineRank[mostSevereCsq ]* dWeig['wRank']
					gpweig+=dWeig['wRank']	
					#~~~~~~~~~~~~~~ find out if a variant is rare 
					listOfRelevantPopulations=[]	# vector with 
					listOfVEPPopulations=["AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_EAS_AF", "gnomAD_FIN_AF", "gnomAD_NFE_AF", "gnomAD_OTH_AF", "gnomAD_SAS_AF"] #data from VEP data format is a float or NULL (empty )
					for vp in listOfVEPPopulations:
						listOfRelevantPopulations.append(dCsq[vp])    

					#listOfAnnovarPopulations=["ExAC_ALL","ExAC_AFR","ExAC_AMR","ExAC_EAS","ExAC_FIN","ExAC_NFE","ExAC_OTH","ExAC_SAS"]
					#for ap in listOfAnnovarPopulations: #annovar mette "." quanod e' vuoto
					#	if dInfo[ap] !=  ".": listOfRelevantPopulations.append(dInfo[ap])			

					freqlist = [float(x) for x in listOfRelevantPopulations if x ] 
				
					rare = gp.checkFreq (freqlist, args.r) # check if it is a rare variant (af<thresh)  
					myres.append(rare)
					if rare== True: gpScore+=1*dWeig['wRare'] ; gpweig+=dWeig['wRare']
					elif rare=="NOB": gpScore+=1*dWeig['wNOB']; gpweig+=dWeig['wNOB']
					
					#~~~~~~~~~~ check if row have Embryo,CellCycle,DDD,GmomAD genes
					embryo = DDD = cellcycle = gnomAD = miscarriages = lethal = essential = False  
				
					###~~~  Gene Ontology  embryo development GO:0009790	
					if re.search("ANN_1", decodedLine): embryo=True; gpScore+=1*dWeig['wEmbryoDev']; gpweig+=dWeig['wEmbryoDev']
					myres.append(embryo)

					###~~~ gene in Loss of Function list 	source? 
					if re.search("ANN_4", decodedLine): gnomAD=True; gpScore+=1*dWeig['wGnomeAD']; gpweig+=dWeig['wGnomeAD']
					myres.append(gnomAD)

					###~~~ gene ontology cell cycle  GO:XXXXX
					if re.search("ANN_3", decodedLine): cellcycle=True; gpScore+=1*dWeig['wCellCycle']; gpweig+=dWeig['wCellCycle']
					myres.append(cellcycle)

					###~~~ DDD
					if re.search("ANN_2", decodedLine): DDD=True; gpScore+=1*dWeig['wDDD'];  gpweig+=dWeig['wDDD']
					myres.append(DDD)
					
					###~~~ List of miscarriages Madhuri 
					if re.search("ANN_5", decodedLine): miscarriages=True; gpScore+=1*dWeig['wMisc']; gpweig+=dWeig['wMisc']
					myres.append(miscarriages)

					###~~~ List of candidate lethal genes
					if re.search("ANN_6", decodedLine): lethal=True; gpScore+=1*dWeig['wLethal']; gpweig+=dWeig['wLethal']
					myres.append(lethal)

                                        ###~~~ List of essential genes
					
					if re.search("essential", decodedLine): essential=True; gpScore+=dWeig['wEssential']*int(dInfo['essential']); gpweig+=dWeig['wEssential']
					myres.append(essential)

					##~~~~~~~~~~ ANNOVAR 
					##~~~ CONSERVATION SCORE 
					# add to myres phyloP7way_vertebrate 
					# addd to header 	
					



  					#~~~~~~~~ FILTER ~~~~~~~~~ before write check parameters	
				
					#if dSOTermFineRank[mostSevereCsq ] > args.i and rare==True and embryo==True or dSOTermFineRank[mostSevereCsq ] > args.i and rare==True and cellcycle==True: #or dSOTermFineRank[mostSevereCsq ] > args.i and rare==True and DDD==True or dSOTermFineRank[mostSevereCsq ] > args.i and rare==True and gnomAD==True:	
					#if dSOTermFineRank[mostSevereCsq ] > args.i and csqAllCount>args.c: #and rare==True and csqAllCount>args.c: 
					myres.append(gpScore) 
					myres.append(gpScore/float(gpweig) )
					filemyres.write("\t".join( map(str, myres ) ) )
					filemyres.write('\n')

			

	fileToWrite=open(args.e, 'w')
	for i in listOfErrors: fileToWrite.write( i )
 

if __name__ == "__main__":
	main()
