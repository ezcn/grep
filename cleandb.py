import csv, sys    
from itertools import izip 

myinputcsvfile=sys.argv[1]

# read the first line as  dictionary keys 
headfile = open(myinputcsvfile, 'r')
firstLine =headfile.readline().split() 
headfile.close() 

# read the file and transpose into lists; every column in a list  
mydirtydata = izip(*csv.reader(open(myinputcsvfile, "rb"), delimiter='\t'))

# 
d_data={}; count= 0 
for column in mydirtydata: 
	d_data[firstLine[count]]=set(column[1:]) 
	count +=1 

for i in d_data : 
	res=[i]
	res+=d_data[i]
	print "\t".join(map(str, res ))
