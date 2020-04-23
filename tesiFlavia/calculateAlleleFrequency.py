from collections import Counter

path_node_id = {}
step_node_id = {}

path_id_ref = ''
with open("/home/flavia/Desktop/MstoGfa/seq_gfa","r") as f:
    for line in f:
        line_list = line.strip('\t').split('\t') 
        if line.startswith('P'):
            path_node_id[line_list[1]] = [x.strip('+') for x in line_list[2].strip('\n').split(',')]   #pathid and stepid
            
            if 'ref' in line_list[1].lower():
                path_id_ref = line_list[1]
        elif line.startswith('S'):
            step_node_id[line_list[1]] = line_list[2].strip('\n')    #id and seq  

            
    num_haplotypes = len(path_node_id) - 1            #I don't count the reference haplotype
    
    for pos in range(len(path_node_id[path_id_ref])):
        tmp_seq_in_pos_list = []
        
        for path_id, step_id_list in path_node_id.items():
            if path_id != path_id_ref:
                #print(path_id, seq_current_step)
                seq_current_step = step_node_id[step_id_list[pos]]
                tmp_seq_in_pos_list.append(seq_current_step)
                
        #print('Pos', pos, '- sequences in this pos', tmp_seq_in_pos_list)
        ATCG_counts_dict = Counter(tmp_seq_in_pos_list)
        #print(ATCG_counts_dict)

        #check the reference base and calculate for each allele the frequency. 
	#In this example there are only biallelic alleles (or reference or a different nucleotide)
	
	#for nucleotide, count in ATCG_counts_dict.items():
            #print(nucleotide, count / num_haplotypes)
	
	
        for nucleotide, count in ATCG_counts_dict.items():
        	if count/num_haplotypes != 1:                        #print only AF !=1
        		print(pos,nucleotide,count/num_haplotypes)


	
		
				
