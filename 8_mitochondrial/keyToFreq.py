mt = open('/home/giuliana/Desktop/Lab/frequency_mt/frequenciesMT.frq')
#mt= open('/home/giuliana/Desktop/Lab/frequency_mt/prova.frq') 
with mt as f:
    next(f)
    for line in f:
         # line split 
         mysplit=line.split()
         #print(line) 
         # generate  ref key 
         mykey_ref=mysplit[0]+":"+ mysplit[1]+":"+ "/" + mysplit[4][0]
     
         # reference data 
         myref = mysplit[4]
         myreffreq =myref.split(":")[1]

         # print mykey, my ref data 
         print(mykey_ref,myreffreq)

         #generate alternate key
         mykey_alt=mysplit[0]+":"+ mysplit[1]+":"+ "/" + mysplit[5][0]
       
         # alternate  data
         myalt = mysplit[5] 
         myaltfreq = myalt.split(":")[1]
        
         #print mykey alternate data 
         print(mykey_alt,myaltfreq)
