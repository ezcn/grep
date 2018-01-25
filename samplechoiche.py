#######
# 



## read fields and subfileds of the data base 

db_fields={}
 
for line in open ("databasefileds.txt"): 
	x=line.rstrip().split()
	field=x[1]; subfield=x[2]
	if not field in db_fields: db_fields[field]=[]
	db_fields[field].append(subfield)  	

#print dbfields

## for relevant fields read tables of scores ; create the dictionary of the scores  

listofrelevantfields=['Info_aborto']
db_score={}
for fie in listofrelevantfields: 
	for sub in db_fields[fie]: 
		filename="tabs/tab_%s" %(sub ) 
		#print filename

		for line in open (filename): 
			y=line.rstrip().split()	
			category=y[0]; score= int(y[1]) 
			db_score[(fie,sub, category)]=score

#print db_score 

######################################
# read sample data and output score 
countlines=0
db_individuals={}
for pline in open ("test.patients"):
	z=pline.split()   
	#print '################'
	if countlines==0 : 
		listoffields=z ; countlines+=1 
	elif countlines==1 : 
		listofsubfields=z ; countlines+=1 
		myreference=zip(listoffields, listofsubfields) 
		#print myreference 
	else: 
		counter=0; db_individuals[countlines]={}
		for temp_category in z : 

			temp_field = myreference[counter][0]
			temp_subfield= myreference[counter][1]

			if not temp_field in db_individuals[countlines]: db_individuals[countlines][temp_field]={}
			
			if  temp_field=="Clinica_anagrafica" :
				db_individuals[countlines][temp_field][temp_subfield]=temp_category 	
			
			elif temp_field in listofrelevantfields: 
				#print (temp_field, temp_subfield, temp_category) 
				db_individuals[countlines][temp_field][temp_subfield]=db_score[(temp_field, temp_subfield, temp_category)]
								
			counter+=1 
		countlines+=1 

##### 
#print db_individuals 
title=[]
anagraphic = db_fields["Clinica_anagrafica"] 
title+=anagraphic
title.append("tot_score")
title+=listofrelevantfields
print "\t".join(title)
 
for indv in db_individuals: 
	res=[]; temp_score=[]
	for dd in db_fields["Clinica_anagrafica"] : res.append(db_individuals[indv]["Clinica_anagrafica"][dd])

	for ff in listofrelevantfields: 
		myscore=sum(db_individuals[indv][ff].values()) 
 		temp_score.append(myscore )

	res.append(sum(temp_score))
	res+=temp_score 
	print "\t".join(map(str, res))  



