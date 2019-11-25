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
	#parser.add_argument("-m", help="path to metadata file",type=str,required=True)
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
##########~~~~~~~~~~~~~~ Read metadata and randomize sample to choose
#	Region=[]
#	Sample=[]
#	
#	for line in open(args.m, 'r'):
#		if re.match('sample', line):
#			header= line.rstrip().split()
#		else:
#			other=line.rstrip().split()
#			dMeta= dict(zip(header, other))
#
#			Region.append(dMeta['region'])
#			Sample.append(dMeta['sample'])
#	
#	dSampleRegion=dict(zip(Sample, Region))
#	
#	EUR = [key  for (key, value) in dSampleRegion.items() if value == 'EUROPE']
#	random.seed(899)
#	
#	sampleToConsider=random.sample(EUR, 6)


##########~~~~~~~~~~~~~~  Loop of vcf lines 
	sys.stdout=open(args.o, 'w')  
	listOfErrors=[]
	dInfo={}
	print('Chr','\t','SNV_mean','\t','Indel_mean','\t','SequenceAlt_mean','\t','Insertion_mean','\t','transcript_ablation','\t','splice_acceptor_variant','\t','splice_donor_variant','\t','stop_gained','\t','frameshift_variant','\t','stop_lost','\t','start_lost','\t','transcript_amplification','\t','inframe_insertion','\t','inframe_deletion','\t','missense_variant','\t','protein_altering_variant','\t','splice_region_variant','\t','incomplete_terminal_codon_variant','\t','start_retained_variant','\t','stop_retained_variant','\t','synonymous_variant','\t','coding_sequence_variant','\t','mature_miRNA_variant','\t','5_prime_UTR_variant','\t','3_prime_UTR_variant','\t','non_coding_transcript_exon_variant','\t','intron_variant','\t','NMD_transcript_variant','\t','non_coding_transcript_variant','\t','upstream_gene_variant','\t','downstream_gene_variant','\t','TFBS_ablation','\t','TFBS_amplification','\t','TF_binding_site_variant','\t','regulatory_region_ablation','\t','regulatory_region_amplification','\t','feature_elongation','\t','regulatory_region_variant','\t','feature_truncation','\t','intergenic_variant')
	element=['SNV','indel','sequence_alteration','deletion','insertion']
	counter=[]
	for i in element:
		counter.append(0)
	counterSOTerm=[]
	for i in lSOTerm:
		counterSOTerm.append(0)
	dcount=dict(zip(element,counter))
	dmean=dict(zip(element,counter))
	dfreq=dict(zip(element,counter))
	dSOTcount=dict(zip(lSOTerm,counterSOTerm))
	dSOTmean=dict(zip(lSOTerm,counterSOTerm))
	dSOTfreq=dict(zip(lSOTerm,counterSOTerm))
	for line in gzip.open(args.f, 'r'):
		decodedLine=line.decode()  ## line.decode() is necessary to read encoded data using gzip in python3
		if re.match('#', decodedLine):
			if re.search("ID=CSQ" , decodedLine ):
				csqHeader=decodedLine.rstrip().split(":")[1].lstrip().rstrip("\">").split("|")
		else:
			linesplit=decodedLine.rstrip().split()
			mychr=linesplit[0]; mypos=linesplit[1]; myref=linesplit[3]; myalt=linesplit[4] ## basic info

			##~~ split INFO field
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
				
				for c in dCsq['Consequence'].split("&"):
					#~~~~~~~~~~~~ assign severity score at the  most severe csq
					myindexes=[]
					myindexes.append(lSOTerm.index(c ))
					mostSevereCsq=lSOTerm[min(myindexes)]
					linesplit[7]=tempinfo	# reset info field 
				for key in dSOTcount:
					if dCsq['Consequence'] in key and freqCSQ_REF_ALT[0]!='NA':
						dSOTcount[dCsq['Consequence']]+=int(1)
						dSOTfreq[dCsq['Consequence']]+=float(freqCSQ_REF_ALT[0])
                                                        #print(dSOTcount,'\t',dSOTfreq)
				for key in dcount:
					if dCsq['VARIANT_CLASS'] in key and freqCSQ_REF_ALT[0]!='NA':
						dcount[dCsq['VARIANT_CLASS']]+=int(1)
						dfreq[dCsq['VARIANT_CLASS']]+=float(freqCSQ_REF_ALT[0])

	for key in dcount:
			if dcount[key]!=0:
				dmean[key]=(dfreq[key]/dcount[key])
	for key in dSOTcount:
			if dSOTcount[key]!=0:
				dSOTmean[key]=(dSOTfreq[key]/dSOTcount[key])


	print(mychr,'\t',dmean['SNV'],'\t',dmean['indel'],'\t',dmean['sequence_alteration'],'\t',dmean['insertion'],'\t',dSOTmean['transcript_ablation'],'\t',dSOTmean['splice_acceptor_variant'],'\t',dSOTmean['splice_donor_variant'],'\t',dSOTmean['stop_gained'],'\t',dSOTmean['frameshift_variant'],'\t',dSOTmean['stop_lost'],'\t',dSOTmean['start_lost'],'\t',dSOTmean['transcript_amplification'],'\t',dSOTmean['inframe_insertion'],'\t',dSOTmean['inframe_deletion'],'\t',dSOTmean['missense_variant'],'\t',dSOTmean['protein_altering_variant'],'\t',dSOTmean['splice_region_variant'],'\t',dSOTmean['incomplete_terminal_codon_variant'],'\t',dSOTmean['start_retained_variant'],'\t',dSOTmean['stop_retained_variant'],'\t',dSOTmean['synonymous_variant'],'\t',dSOTmean['coding_sequence_variant'],'\t',dSOTmean['mature_miRNA_variant'],'\t',dSOTmean['5_prime_UTR_variant'],'\t',dSOTmean['3_prime_UTR_variant'],'\t',dSOTmean['non_coding_transcript_exon_variant'],'\t',dSOTmean['intron_variant'],'\t',dSOTmean['NMD_transcript_variant'],'\t',dSOTmean['non_coding_transcript_variant'],'\t',dSOTmean['upstream_gene_variant'],'\t',dSOTmean['downstream_gene_variant'],'\t',dSOTmean['TFBS_ablation'],'\t',dSOTmean['TFBS_amplification'],'\t',dSOTmean['TF_binding_site_variant'],'\t',dSOTmean['regulatory_region_ablation'],'\t',dSOTmean['regulatory_region_amplification'],'\t',dSOTmean['feature_elongation'],'\t',dSOTmean['regulatory_region_variant'],'\t',dSOTmean['feature_truncation'],'\t',dSOTmean['intergenic_variant'])


if __name__ == "__main__":
	main()
