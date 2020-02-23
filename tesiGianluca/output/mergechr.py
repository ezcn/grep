import glob ,sys 

"""
usage :python3 mergechr.py 'snv1/*' csqimpact.tsv > freq.all.SNV1.out
python3 mergechr.py 'snv05/*' csqimpact.tsv   > freq.all.SNV05.out
"""
pathToFolder=sys.argv[1]
csqimpactFile=sys.argv[2]

def combineMeanSD (listOfReplicates): 
	"""
	https://www.statstodo.com/CombineMeansSDs_Pgm.php
	listOfReplicates=( {'n': 10, 'm': 11.8, 'sd': 2.4},  {'n': 20, 'm': 15.3, 'sd': 3.2},  {'n': 15, 'm': 8.4, 'sd': 4.1})
	"""
	tn=sum([rep['n'] for rep in listOfReplicates])
	tx= sum([rep['m'] * rep['n'] for rep in listOfReplicates ] ) 
	txx= sum([rep['sd']**2 *(rep['n']-1) + (rep['m'] * rep['n'])**2/rep ['n'] for rep in listOfReplicates]) 
	combMean=tx/tn
	combSD=((txx-tx**2/tn) /(tn-1))  **(1/2)
	return combMean, combSD 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
soD={}
for line in open (csqimpactFile): 
    y=line.rstrip().split('\t') 
    soD[y[0]]=(y[3],y[4])  # soD[type]=(Display term, IMPACT) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
merged={}
for ty in soD: merged[ty]=[]
path = pathToFolder
files=glob.glob(path)   
#print(files) 

for f in files: 
    firstLine=True 
    for line in open (f) : 
        x=line.rstrip().split()
        if firstLine: 
            header=x; firstLine=False 
        else: 
            tempD=dict(zip(header, x ))
            if int(tempD['n'])>5:
                typeD={'n': float(tempD['n']), 'm': float(tempD['muZFreq']), 'sd': float(tempD['sdZFreq'])}
                merged[tempD['type']].append(typeD) 

print ("\t".join(['type', 'display', 'impact', 'n', 'Zmu', 'Zsd'])) 
for typ in merged:
    if len(merged[typ])> 0 : 
        #print ('###################################')
        totalN=sum([n['n'] for n in merged[typ] ])
        display =soD[typ][0]; impact=soD[typ][1]
        res=combineMeanSD( merged[typ])
        print('\t'.join(map(str , [typ, display,impact ,  totalN, res[0], res[1]])) )
        #print (len(merged[typ]) , merged[typ])
