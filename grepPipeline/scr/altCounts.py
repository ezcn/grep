import re, sys, argparse
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", help="path to input file  ",type=str, required= True)
	parser.add_argument("-id", help="sample id  ",type=str, required= True)
	args = parser.parse_args()
	
	for line in open(args.i, 'r'):
		if re.match('CHROM', line):
			print('key' + '\t' + 'ID' + '\t' + 'ALTcount')
		else:
			a=line.rstrip().split()
			new=a[5].split(":")
			print(a[0]+ ':' + a[1] + ':/' + new[0] + '\t' + args.id + '\t' + new[1] )
			 
if __name__ == "__main__":
	main()		  
