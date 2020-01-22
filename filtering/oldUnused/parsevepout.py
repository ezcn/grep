## Take as input the output of VEP and outputs a table with chr pos allele *score*
## *score* is based on VEP consequence rank need to be IMproved 


myinputvep=open("testvep.tsv", "r") ## tab delimited out of vepcsq
myvepcsq=open ("csqimpact.tsv", "r") ## https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html

##  1) READ VEP consequences rank ########
## read VEP rank and create dSOTermRank: dictionary with score assigend to SOTerm based on dRank  *IMPROVE scoring*  
dRank={"HIGH":4, "LOW": 2, "MODERATE":3, "MODIFIER":1}
dSOTermRank={} 

countlinesCsq= True
for csqLine in myvepcsq:
        if countlinesCsq:
                csqTitle=csqLine.rstrip().split('\t')
                countlinesCsq=False
        else:
                myCsqList=csqLine.rstrip().split('\t')
                dCsq= dict(zip(csqTitle, myCsqList ))  ## create a temporary dictionary with the title column and the content of the row 
                #print("##########################" )
                #print (dCsq)
                dSOTermRank[dCsq['SO term']]=dRank[dCsq['IMPACT']]
                #print('~~~~~~~~~~~~~~~~~')
#print( dSOTermRank)

## 2) Process VEP output
countlines=0
for line in myinputvep:
        if countlines ==0 :
                title =line.rstrip().split()
                countlines+=1
        else:
                myrow=line.rstrip().split()
                drow= dict(zip(title, myrow ))   ## create a temporary dictionary with the title column and the content of the row 
                #print(drow)
                print ( drow['Chr'], drow['Pos'], drow['Allele'], dSOTermRank[drow['Consequence'] ])
