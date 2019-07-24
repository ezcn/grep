import re 
import greplib as gp 

myvcfinput=open("AS064.test.vcf") 

########################################################


myvepcsq=open ("csqimpact.tsv", "r")

##  READ VEP consequences rank ########
dRank={"HIGH":4, "LOW": 2, "MODERATE":3, "MODIFIER":1}
dSOTermRank={}
lSOTerm=[]

countlinesCsq= True
for csqLine in myvepcsq:
	if countlinesCsq:
		csqTitle=csqLine.rstrip().split('\t')
		countlinesCsq=False
	else:
		myRowList=csqLine.rstrip().split('\t')
		dCsq= dict(zip(csqTitle, myRowList ))
		dSOTermRank[dCsq['SO term']]=dRank[dCsq['IMPACT']]
		lSOTerm.append(myRowList[0])

#print (lSOTerm)


#########################################################

dInfo={}
header=["chr", "pos", "csqAllel", "csqAlleleCount", "GTLiklihood" , "ENSTID", "ImpactScore", "rare","Embryo","GnomAD","CellCycle","DDD"]
print("\t".join(map(str, header) )  ) 

for line in myvcfinput: 
	if not re.match("#", line): 
		#print("this is a new line ") 
		## line split by  tab 
		linesplit=line.rstrip().split()
		
		## basic info 
		mychr=linesplit[0]; mypos=linesplit[1]; myref=linesplit[3]; myalt=linesplit[4] 

		## split INFO field
		tempinfo=linesplit[7] 
		for i in tempinfo.split(";"):  
			temp=i.split("=") 
			dInfo[temp[0]]=temp[1]

		## split FORMAT field
		tempformattitle=linesplit[8].split(":")
		tempformatcontent=linesplit[9].split(":")
		dFormat=dict(zip(tempformattitle, tempformatcontent))

		## work on dInfo[CSQ]
		## split for multiple consequences separated by ","
		multipleCsq=dInfo["CSQ"].split(",") 

		for mcsq in multipleCsq:  ### single consequence  
			myres=[]
			myres+=[mychr, mypos]
			dCsq=dict(zip(csqHeader, mcsq.split("|") ))  #############    ALL VEP INFO 
			
			#~~~~~~~~~
			mycsqAllele=dCsq["Allele"] ## identify the allele with consequences 

			#~~~~~~~~~~~
			myres+= gp.csqAlleleFeatures(mycsqAllele, myalt, int(dInfo["AC"]), dFormat["GL"] ) ## features of csqAll
		
			#~~~~~~~~~~~~~
			myres.append(dCsq['Feature'])
				
			#~~~~~~~~~~~~~~~~			
			myind=[]
			for tl in dCsq['Consequence'].split("&"): 
				myind.append(lSOTerm.index(tl ))	
			mostSevereCsq=lSOTerm[min(myind)]
			#print(dCsq['Consequence']) 
			#print(mostSevereCsq)
			myres.append( dSOTermRank[mostSevereCsq ]) ## score based on the impact of the consequence       	

			#~~~~~~~~~~~~~~~~~~~~~
			thresh=0.01
			freqlist = [float(x) for x in  [dCsq["AFR_AF"],dCsq["AMR_AF"],dCsq["EAS_AF"],dCsq["EUR_AF"],dCsq["SAS_AF"]] if x ] 
			rare = gp.checkFreq (freqlist, thresh) # check if it is a rare variant (af<thresh) 
			myres.append(rare)
			#print ( "\t".join( map(str, myres) )  )					
			#~~~~~~~~~~~~~~~~~~~~~~~~~~
			
			if re.search("annotation", line): myres.append("Embryo")
			else: myres.append("na")
			if re.search("ANN3", line): myres.append("GnomAD")
			else: myres.append("na")
			if re.search("ANN2", line): myres.append("CellCycle")
			else: myres.append("na")
			if re.search("ANN1", line): myres.append("DDD")
			else: myres.append("na")
							
							
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
			print ( "\t".join( map(str, myres) )  )

	else: 
		if re.search("ID=CSQ" ,line ): 
			csqHeader=line.rstrip().split(":")[1].lstrip().rstrip("\">").split("|")		
			#print (csqHeader)	



 



