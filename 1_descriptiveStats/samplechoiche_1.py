####### 

#dbfields_file="/home/enza/oogaprotocol/IMMA/abortion_db/AbortionAS-AVdatabase-IMMA05-02-2018.dbfileds.tsv"
dbfields_file="/home/enza/oogaprotocol/IMMA/abortion_db/AbortionAS-AVdatabase-IMMA.09-02-2018.dbfields.tsv"

#patient_data_file="/home/enza/oogaprotocol/IMMA/abortion_db/AbortionAS-AVdatabase-IMMA05-02-2018.claning.tsv"
patient_data_file="/home/enza/oogaprotocol/IMMA/abortion_db/AbortionAS-AVdatabase-IMMA.09-02-2018.minusfirstline.tsv"


listofIDfields=['Type', 'ID', 'Nascita', 'CodiceFiscale', 'Interview_date']
listofrelevantsub=['Menarche_age', 'BMI', 'Miscarriage']
db_score={}
#for fie in listofrelevantfields: 
for sub in listofrelevantsub: 
	subfield_filename="tabs/tab_%s" %(sub ) 
	#print filename

	with open(subfield_filename) as f:
    		next(f)
 		for line in f:
			y=line.rstrip().split()	
			category=y[0]; score= int(y[1]) 
			db_score[(sub, category)]=score

#print db_score 

######################################
# read sample data and output score 
countlines=0
db_individuals={}
for pline in open (patient_data_file):
	z=pline.split('\t')  
 
	if countlines==0 : 
		listofsubfields=z ; countlines+=1 
		#print listofsubfields 
	else: 
		countlines+=1
		count_values=0; db_individuals[countlines]={}
		for temp_value in z :  # iterates over individual's values 
			temp_subfield= listofsubfields[count_values] ; count_values+=1
			#print temp_subfield	
			if  temp_subfield in listofIDfields:
				db_individuals[countlines][temp_subfield]=temp_value 	
			
			elif temp_subfield in listofrelevantsub: 
				db_individuals[countlines][temp_subfield]=db_score[(temp_subfield, temp_value)]
								

##### 
#print db_individuals 
title=[]
title+=listofIDfields
title.append("tot_score")
title+=listofrelevantsub
print "\t".join(title)
 
for indv in db_individuals: 
	res=[]; temp_score=[]
	for dd in listofIDfields: res.append(db_individuals[indv][dd])

	for ff in listofrelevantsub: 
		myscore=db_individuals[indv][ff]
 		temp_score.append(myscore )

	res.append(sum(temp_score))
	res+=temp_score 
	print "\t".join(map(str, res))  



