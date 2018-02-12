
myinfile="/home/enza/oogaprotocol/IMMA/abortion_db/AbortionAS-AVdatabase-IMMA05-02-2018.cleaning.minusfirstline.tsv"

typesofabortion=[]
countline=0 
d_data={}
for line in open (myinfile) : 
	x=line.split("\t")  
	if countline > 0 : 
		mytype=x[0]; myID=x[1]; outcome=[]
		for i in range(74, 160, 11): 
			if not x[i] in typesofabortion: typesofabortion.append(x[i]) 
			outcome.append(x[i]) 
		d_data[(mytype, myID)]=outcome 
	else: countline +=1

title=["nb_type", "nb_ID"]
for ttt in typesofabortion: 
	if not ttt=="na": title.append(ttt) 
print "\t".join(title)


for myid in d_data: 
	res=list(myid)
	for tt in typesofabortion: 
		if not tt=="na": 
			res.append(d_data[myid].count(tt))
	print "\t".join(map (str, res) )	
