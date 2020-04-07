import statistics
import matplotlib.pyplot as plt
 

mean_fst_list = []
for num_rep in range(1,101):   #num replicati
    f1 = open('/home/flavia/Desktop/ms_to_vcf/simul_reali/100Replicate/frqAlPop1/ms_rep{}.bcf.pop1.vcf.frq'.format(num_rep)) #Fst all0 e all1 pop1
    f2 = open('/home/flavia/Desktop/ms_to_vcf/simul_reali/100Replicate/frqAlPop2/ms_rep{}.bcf.pop2.vcf.frq'.format(num_rep)) #Fst all0 e all1 pop2

    fst_list = []
    f1.readline() # Skip header
    f2.readline() # Skip header
    for (line1, line2) in zip(f1, f2):
        freq0_pop1 = line1.strip('\n').split('\t')[-2]   #penultima colonna
        freq0_pop2 = line2_list = line2.strip('\n').split('\t')[-2]  #penultima colonna
        freq0_pop1 = float(freq0_pop1.split('0:')[1])
        freq0_pop2 = float(freq0_pop2.split('0:')[1])

        mean = statistics.mean([freq0_pop1, freq0_pop2])
        #print(mean)
        fst = statistics.pvariance([freq0_pop1, freq0_pop2])/((mean)*(1-mean))
    #except ZeroDivisionError:
        #z = 0
        #print(freq0_pop1, freq0_pop2)
        fst_list.append(fst)
       	
    f1.close()
    f2.close()


    mean_fst_list.append(statistics.mean(fst_list))

print(mean_fst_list)

#plt.style.use('ggplot')
#plt.hist(mean_fst_list, bins=10)
#plt.show()