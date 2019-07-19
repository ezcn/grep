myinputvep=open("testvep.tsv", "r")
myvepcsq=open ("csqimpact.tsv", "r")

##  READ VEP consequences rank ########
dRank={"HIGH":4, "LOW": 2, "MODERATE":3, "MODIFIER":1}
dSOTermRank={}

countlinesCsq= True
for csqLine in myvepcsq:
        if countlinesCsq:
                csqTitle=csqLine.rstrip().split('\t')
                countlinesCsq=False
        else:
                myCsqList=csqLine.rstrip().split('\t')
                dCsq= dict(zip(csqTitle, myCsqList ))
                #print("##########################" )
                #print (dCsq)
                dSOTermRank[dCsq['SO term']]=dRank[dCsq['IMPACT']]
                #print('~~~~~~~~~~~~~~~~~')
#print( dSOTermRank)

###   Process VEP output
countlines=0
for line in myinputvep:
        if countlines ==0 :
                title =line.rstrip().split()
                countlines+=1
        else:
                myrow=line.rstrip().split()
                drow= dict(zip(title, myrow ))
                #print(drow)
                print ( drow['Chr'], drow['Pos'], drow['Allele'], dSOTermRank[drow['Consequence'] ])
