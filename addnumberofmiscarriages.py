#!/usr/bin/python

infile=open("/home/enza/oogaprotocol/IMMA/abortion_db/AbortionAS-AVdatabase-IMMA05-02-2018.nbofmiscarriages.tsv") 

dbfile="/home/enza/oogaprotocol/IMMA/abortion_db/AbortionAS-AVdatabase-IMMA05-02-2018.cleaning.tsv"

#line = infile.readline()
d_nb= dict()
countline=0
for  line in infile: 
	content = line.rstrip().split()
	if countline==0: 
		countline+=1 
		title=content
	else: 
		myid = content[1]
		d_nb[myid]=content

countdline=0 
for dline in open(dbfile): 
	dline = dline.rstrip().split("\t")  
	dbid=dline[1]	
	if countdline==0:
		countdline+=1 
		for i in range(len(title)) :  
			dline.append("Previous_summary")
		#print "\t".join(dline )
	elif countdline==1 : 
		countdline+=1  
		dline+=title 
		#print "\t".join(dline )
	else: 
		dline+=d_nb[dbid]	

	print "\t".join(dline )
