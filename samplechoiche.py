####### 

#dbfields_file="/home/enza/oogaprotocol/IMMA/abortion_db/AbortionAS-AVdatabase-IMMA05-02-2018.dbfileds.tsv"
dbfields_file="/home/enza/oogaprotocol/IMMA/abortion_db/AbortionAS-AVdatabase-IMMA.09-02-2018.dbfields.tsv"

#patient_data_file="/home/enza/oogaprotocol/IMMA/abortion_db/AbortionAS-AVdatabase-IMMA05-02-2018.claning.tsv"
patient_data_file="/home/enza/oogaprotocol/IMMA/abortion_db/AbortionAS-AVdatabase-IMMA.09-02-2018.tsv"



## read fields and subfileds of the data base 

db_fields={}
 
for line in open (dbfields_file): 
	x=line.rstrip().split()
	field=x[1]; subfield=x[2]
	if not field in db_fields: db_fields[field]=[]
	db_fields[field].append(subfield)  	

#print dbfields

## for relevant fields read tables of scores ; create the dictionary of the scores  

listofrelevantfields=['Patient']
listofrelevantsub=['Menarche_age', 'BMI']
db_score={}
for fie in listofrelevantfields: 
	for sub in db_fields[fie]: 
		if sub in listofrelevantsub: 
			subfield_filename="tabs/tab_%s" %(sub ) 
			#print filename

			with open(subfield_filename) as f:
    				next(f)
 				for line in f:
			#for line in open (filename): 
					y=line.rstrip().split()	
					category=y[0]; score= int(y[1]) 
					db_score[(fie,sub, category)]=score

#print db_score 

######################################
# read sample data and output score 
countlines=0
db_individuals={}
for pline in open (patient_data_file):
	z=pline.split()   
	#print '################'
	if countlines==0 : 
		listoffields=z ; countlines+=1 
	elif countlines==1 : 
		listofsubfields=z ; countlines+=1 
		myreference=zip(listoffields, listofsubfields) 
		print myreference 
	else: 
		counter=0; db_individuals[countlines]={}
		for temp_value in z :  # iterates over individual's values 

			temp_field = myreference[counter][0]
			temp_subfield= myreference[counter][1]

			#if not temp_field in db_individuals[countlines]: db_individuals[countlines][temp_field]={}
			
			if  temp_field=="Patient_ID" :
				db_individuals[countlines][temp_field][temp_subfield]=temp_value 	
			
			elif temp_field in listofrelevantfields: 
				if not temp_field in db_individuals[countlines]: db_individuals[countlines][temp_field]={}
				if temp_subfield in listofrelevantsub: 
				#print (temp_field, temp_subfield, temp_category) 
					db_individuals[countlines][temp_field][temp_subfield]=db_score[(temp_field, temp_subfield, temp_value)]
								
			counter+=1 
		countlines+=1 

##### 
#print db_individuals 
title=[]
anagraphic = db_fields["Patient_ID"] 
title+=anagraphic
title.append("tot_score")
title+=listofrelevantfields
print "\t".join(title)
 
for indv in db_individuals: 
	res=[]; temp_score=[]
	for dd in db_fields["Patient_ID"] : res.append(db_individuals[indv]["Patient_ID"][dd])

	for ff in listofrelevantfields: 
		myscore=sum(db_individuals[indv][ff].values()) 
 		temp_score.append(myscore )

	res.append(sum(temp_score))
	res+=temp_score 
	print "\t".join(map(str, res))  



