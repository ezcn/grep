import csv, sys, os     
from itertools import izip 
from pathlib import Path


#myinputcsvfile=sys.argv[1]
#myinputcsvfile="/home/enza/oogaprotocol/IMMA/abortion_db/AbortionAS-AVdatabase-IMMA05-02-2018.cleaning.minusfirstline.tsv"
myinputcsvfile="/home/enza/oogaprotocol/IMMA/abortion_db/AbortionAS-AVdatabase-IMMA.09-02-2018.minusfirstline.tsv"


#mydbfields="/home/enza/oogaprotocol/IMMA/abortion_db/AbortionAS-AVdatabase-IMMA05-02-2018.dbfileds.tsv"
mydbfields="/home/enza/oogaprotocol/IMMA/abortion_db/AbortionAS-AVdatabase-IMMA.09-02-2018.dbfields.tsv"


# read the first line as  dictionary keys 
headfile = open(myinputcsvfile, 'r')
firstLine =headfile.readline().rstrip().split("\t") 
headfile.close() 

# read the file and transpose into lists; every column in a list  
mydirtydata = izip(*csv.reader(open(myinputcsvfile, "rb"), delimiter='\t'))

#for mm in mydirtydata: print mm 
# 
d_data={}; count= 0 
for column in mydirtydata: 
	d_data[firstLine[count]]=set(column[1:]) 
	count +=1 


#print d_data['Type_2_diabetes_non-insulin-dependent'] 
####### make the tables 
tabdir=Path("tabs")
 
if not tabdir.is_dir(): tabdir.mkdir()  #os.mkdir(tabdir)

 
for dbfield  in d_data :
	tab_file="%s/tab_%s" %(tabdir, dbfield)
	tab_file_path = Path(tab_file)
	if not tab_file_path.is_file():
		tab_out=open(tab_file, "w")
		sys.stdout=tab_out  
		print dbfield, "score"
		for item in d_data[dbfield]: print item  
		

## print a list
mylistofvalues_file=open("dbfields.content.tsv", "w")
sys.stdout=mylistofvalues_file 
for line in open (mydbfields) :
	dline=line.rstrip().split("\t")
	myfield=dline[2]
	for item in d_data[myfield]: dline.append(item)  	 
	print "\t".join(dline) 



