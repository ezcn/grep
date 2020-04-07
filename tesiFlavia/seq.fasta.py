with open("/home/flavia/Desktop/MstoGfa/seqwa100", "r") as f:
	next(f)
	filedata=f.readlines()
	columns = [1]

	

	seq = []
	
	for line in filedata:
	    columns = line.split()
	    seq.append(columns[1])
	    

value = []
for line in filedata:
	columns = line.split()
	value.append(columns[0])
	    
x = list(zip(value,seq))


ofile = open("seq.fasta100", "w")

for i in range(len(x)):
	ofile.write(">" + value[i] + "\n" +seq[i] + "\n")


ofile.close()





	 
	      


''' 
inp = open('/home/flavia/tools/seqfilewa','r')
linesinp=inp.readlines()
inp_column_number = 1 

resultinp=[] 
for x in linesinp:
    resultinp.append(x.split())
    #resultinp.append(x.rstrip("\n"))
    print("> " +  x + "\n")
inp.close()
#print(resultinp)



inp = ("/home/flavia/tools/seqfilewa", "r")
outp = open('fasta.fa', "w")
print ("Convertion")
	
for a in inp:
	a = a.rstrip("\n")
	outp.write("> " +  1 + "\n")
	#outp.write([11:-4] + "\n")
	#print ("Done")

	#inp.close()
	#outp.close()
'''