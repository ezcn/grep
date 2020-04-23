ref_seq = ''

path_output = ('seqwa.seqgen')

with open(path_output) as f:
    num_ref, len_seq = [int(x) for x in f.readline().strip(' \n').split(' ')] 
    num_ref += 1

    for line in f:
        num, seq = line.strip('\n').split('\t')

        if int(num) == num_ref:
            ref_seq = seq
            break

if ref_seq == '':
    print('No reference found!')
    exit()


# The first step_id is the reference step id
pos_to_step_ids_in_that_pos_dict = {}
    
step_id_to_seq_dict = {}
path_to_nodes_dict = {} # Name and Step List
path_nodes_dict = {}
links_to_write = ''

path_nodes_dict['AncestralReference'] = []
for step_id, nt in enumerate(ref_seq):
    step_id_to_seq_dict[step_id] = nt
    
    path_nodes_dict['AncestralReference'].append(str(step_id))
    
    pos_to_step_ids_in_that_pos_dict[step_id] = [step_id]
    
last_step_id = step_id + 1

with open(path_output) as f:
    f.readline() # Skip line

    stuff_already_seen = set()
    for line in f:
        num, seq = line.strip('\n').split('\t')

        # If it is not the reference...
        if int(num) != num_ref:
            # ...find differences
            
            path_nodes_dict[num] = []
            last_chosen_node = -1
            for pos, (x, y) in enumerate(zip(ref_seq, seq)):
                choosen_step_id = -1
                if x == y:
                    choosen_step_id = pos_to_step_ids_in_that_pos_dict[pos][0] # The first step id is the reference step id
                else:
                    for possible_step_id in pos_to_step_ids_in_that_pos_dict[pos]:
                        if step_id_to_seq_dict[possible_step_id] == y:
                            # Already available a node with the same sequence
                            choosen_step_id = possible_step_id
                            break
                    
                    if choosen_step_id == -1:
                        # I need to create a new step_id
                        step_id_to_seq_dict[last_step_id] = y
                        pos_to_step_ids_in_that_pos_dict[pos].append(last_step_id)
                        choosen_step_id = last_step_id
                    
                        last_step_id += 1
                path_nodes_dict[num].append(str(choosen_step_id))
                
                # Links
                if last_chosen_node != -1:
                    nodo_1 = last_chosen_node
                    nodo_2 = choosen_step_id
                    nodo_12 = str(nodo_1) + '_' + str(nodo_2)
                    if nodo_12 not in stuff_already_seen:
                        stuff_already_seen.add(nodo_12)
                        links_to_write += 'L\t' + str(nodo_1) + '\t+\t' + str(nodo_2) + '\t+\t0M\n'
                        
                last_chosen_node = choosen_step_id
                

                
for step_id, seq in step_id_to_seq_dict.items():
    #print('\t'.join(['S', str(step_id), seq]))
    steps_to_write += 'S\t'+ str(step_id) + '\t' + seq + '\n'
    

for path, node_list in path_nodes_dict.items():
    paths_to_write += ('\t'.join(['P', path,'+,'.join(node_list) + '+', '\n']))

    #for step_id in paths_to_write:
        #print(len(seq))
    

    path_gfa = 'seq_gfa' 
with open(path_gfa, 'w') as fw:
    fw.write(steps_to_write + paths_to_write + links_to_write)
    print(path_gfa + ' written')
