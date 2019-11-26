import re 
import sys
#sys.path.append('/mpba0/vcolonna/gianluca/pythonScript/greplib.py')
sys.path.append('/home/enza/ezcngit/grepl/filtering/greplib.py')
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
	parser.add_argument("-n", help="number of cycles/replicates",type=int,required=True))
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
			dvepCsq= dict(zip(csqTitle, myRowList ))
			dSOTermRank[dvepCsq['SO term']]=dRank[dvepCsq['IMPACT']]
			lSOTerm.append(myRowList[0])

##########~~~~~~~~~~~~~~ Read metadata and randomize sample to choose
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
	
	sampleToConsider=random.sample(EUR, args.n)


##########~~~~~~~~~~~~~~  Loop of vcf lines 
	sys.stdout=open(args.o, 'w')  
	listOfErrors=[]
	preHeader=['Chr','SNV_mean','Indel_mean','SequenceAlt_mean','Insertion_mean']
	Header=preHeader+lSOTerm
	print ("\t".join([i for i in Header]))	
	dInfo={}
	element=['SNV','indel','sequence_alteration','deletion','insertion']
	#counter=[]
	#for i in element:
	#	counter.append(0)
	#counterSOTerm=[]
	#for i in lSOTerm:
	#	counterSOTerm.append(0)
	#dcount=dict(zip(element,counter))
	#dmean=dict(zip(element,counter))
	#dfreq=dict(zip(element,counter))
	#dSOTcount=dict(zip(lSOTerm,counterSOTerm))
	#dSOTmean=dict(zip(lSOTerm,counterSOTerm))
	#dSOTfreq=dict(zip(lSOTerm,counterSOTerm))
	dElem={}; dSOT={}; dImpact={}
	for line in gzip.open(args.f, 'r'):
		decodedLine=line.decode()  ## line.decode() is necessary to read encoded data using gzip in python3
		if re.match('#', decodedLine):
			if re.search("ID=CSQ" , decodedLine ):
				csqHeader=decodedLine.rstrip().split(":")[1].lstrip().rstrip("\">").split("|")
		else:
			linesplit=decodedLine.rstrip().split()
			mychr=linesplit[0]; mypos=linesplit[1]; myref=linesplit[3]; myalt=linesplit[4] ## basic info

			##~~ split INFO field and make a dictionary out of it 
			tempinfo=linesplit[7] 
			for i in tempinfo.split(";"):
				if re.search('=', i):
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
				GTfields=linesplit[9:]
				nbAploidSamples=len(GTfields)*2
				freqCSQ_REF_ALT=gp.AnnotateFreqCSQ_REF_ALT(mycsqAllele,myref, myalt, nbAploidSamples, GTfields)
				
				#~~~~~~~~~~     report for every consequence 
				for c in dCsq['Consequence'].split("&"):
					#~~~~~~~~~~~~ assign severity score at the  most severe csq
					myindexes=[]
					myindexes.append(lSOTerm.index(c ))
					mostSevereCsq=lSOTerm[min(myindexes)]
					linesplit[7]=tempinfo	# reset info field 
				
				if freqCSQ_REF_ALT[0]!='NA': 
	
					if not dCsq['Consequence'] in dSOT: dSOT[dCsq['Consequence']]=[0,0] #inizialize di dictionary with [counter, allele freq] if the key is not present 
					else: 
						dSOT[dCsq['Consequence']][0]+=1 #add +1 to the counter 	
						dSOT[dCsq['Consequence']][1]+=float(freqCSQ_REF_ALT[0]) #add the value of the consequence allele 

					if not dCsq['VARIANT_CLASS'] in dElem: dElem[dCsq['VARIANT_CLASS']]=[0,0]
					else: 
						dElem[dCsq['VARIANT_CLASS']][0]+=1		
						dElem[dCsq['VARIANT_CLASS']][1]+=float(freqCSQ_REF_ALT[0]) 
	
				else: listOfErrors.append((mychr, mypos,myref, myalt, dCsq["Allele"]) ) #to be printed in the error file to compare allele matching  
				

#				for key in dSOTcount:
#					if dCsq['Consequence'] in key and freqCSQ_REF_ALT[0]!='NA':
#						dSOTcount[dCsq['Consequence']]+=int(1)
#						dSOTfreq[dCsq['Consequence']]+=float(freqCSQ_REF_ALT[0])
                                                        #print(dSOTcount,'\t',dSOTfreq)
#				for key in dcount:
#					if dCsq['VARIANT_CLASS'] in key and freqCSQ_REF_ALT[0]!='NA':
#						dcount[dCsq['VARIANT_CLASS']]+=int(1)
#						dfreq[dCsq['VARIANT_CLASS']]+=float(freqCSQ_REF_ALT[0])
#
#	for key in dcount:
#			if dcount[key]!=0:
#				dmean[key]=(dfreq[key]/dcount[key])
#	for key in dSOTcount:
#			if dSOTcount[key]!=0:
#				dSOTmean[key]=(dSOTfreq[key]/dSOTcount[key])
 
	myres=[mychr]
	for vel in element:
		if vel in dElem:  myres.append(dElem[vel][1]/float(dElem[vel][0]))
		else: myres.append('na')  

	for vcsq in lSOTerm: 
		if vcsq in dSOT: myres.append(dSOT[vel][1]/float(dSOT[vel][0]))
		else: myres.append('na')

	print ("\t".join(map(str, myres) ))


#	print(mychr,'\t',dmean['SNV'],'\t',dmean['indel'],'\t',dmean['sequence_alteration'],'\t',dmean['insertion'],'\t',dSOTmean['transcript_ablation'],'\t',dSOTmean['splice_acceptor_variant'],'\t',dSOTmean['splice_donor_variant'],'\t',dSOTmean['stop_gained'],'\t',dSOTmean['frameshift_variant'],'\t',dSOTmean['stop_lost'],'\t',dSOTmean['start_lost'],'\t',dSOTmean['transcript_amplification'],'\t',dSOTmean['inframe_insertion'],'\t',dSOTmean['inframe_deletion'],'\t',dSOTmean['missense_variant'],'\t',dSOTmean['protein_altering_variant'],'\t',dSOTmean['splice_region_variant'],'\t',dSOTmean['incomplete_terminal_codon_variant'],'\t',dSOTmean['start_retained_variant'],'\t',dSOTmean['stop_retained_variant'],'\t',dSOTmean['synonymous_variant'],'\t',dSOTmean['coding_sequence_variant'],'\t',dSOTmean['mature_miRNA_variant'],'\t',dSOTmean['5_prime_UTR_variant'],'\t',dSOTmean['3_prime_UTR_variant'],'\t',dSOTmean['non_coding_transcript_exon_variant'],'\t',dSOTmean['intron_variant'],'\t',dSOTmean['NMD_transcript_variant'],'\t',dSOTmean['non_coding_transcript_variant'],'\t',dSOTmean['upstream_gene_variant'],'\t',dSOTmean['downstream_gene_variant'],'\t',dSOTmean['TFBS_ablation'],'\t',dSOTmean['TFBS_amplification'],'\t',dSOTmean['TF_binding_site_variant'],'\t',dSOTmean['regulatory_region_ablation'],'\t',dSOTmean['regulatory_region_amplification'],'\t',dSOTmean['feature_elongation'],'\t',dSOTmean['regulatory_region_variant'],'\t',dSOTmean['feature_truncation'],'\t',dSOTmean['intergenic_variant'])


if __name__ == "__main__":
	main()
