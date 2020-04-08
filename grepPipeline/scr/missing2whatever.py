#!/usr/bin/python3
import re, sys, argparse, gzip 

'''
usage
python3 missing2whatever.py  -v 1000GP.only-GXYLT1.sitesGrep.vcf.gz  -m . -r '0|0'
REMEMBER TO  WRITE '0|0'   and not 0|0   

'''

def main():
	parser = argparse.ArgumentParser()
	#~~~ input files 
	parser.add_argument("-vcf", help="path to  bgzipped vcf  file ",type=str,required=True)
	parser.add_argument("-r", help="string to replace missing (e.g. 0|0 or 0/0) ",type=str,required=True)
	parser.add_argument("-m", help="format of the missing data in the vcf",type=str,required=True)
	args = parser.parse_args()

	for line in gzip.open(args.vcf): 
		decodedLine=line.decode().rstrip()
		if re.match('#', decodedLine): 
			print(decodedLine) 
		else: 
			x=decodedLine.split()
			y=[ args.r  if re.findall(args.m, i )  else i for i in x[9:] ]
			print ('\t'.join(map(str, x[0:8]+y )))
		
if __name__ == "__main__":
    main()
