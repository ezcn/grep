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
