import os

replicate_to_pos_locus_to_sample_to_genotype_dict = {}

ms_command = ''

with open('/home/flavia/Desktop/MstoGfa/out2pop.ms') as f:
	ms_command = f.readline().strip()

	# Skit the following 1 row(s)
	for _ in range(1):
		f.readline()

	# Fast implementation (to improve for huge simulations)
	for num_rep, replicate in enumerate(f.read().split('//')[1:]):
		print('Replicate num.', num_rep)

		replicate_list = replicate.strip().split('\n')

		# Params
		param_dict = {}
		for num_row, row in enumerate(replicate_list):
			if ':' not in row:
				break
			param_name, param_info = row.split(': ')
			if param_name in row:
				param_dict[param_name] = param_info


		if 'positions' not in param_dict.keys():
			print('ERROR: no positions found for the replicate num. {}'.format(num_rep))
			continue
		param_dict['positions'] = param_dict['positions'].strip().split(' ')

		if 'segsites' not in param_dict.keys():
			print('ERROR: no segsites found for the replicate num. {}'.format(num_rep))
			continue
		else:
			param_dict['segsites'] = int(param_dict['segsites'])

		if 'prob' not in param_dict.keys():
			print('WARNING: no prob found for the replicate num. {}'.format(num_rep))
		else:
			param_dict['prob'] = float(param_dict['prob'])

		

		steps_to_write = ''
		paths_to_write = ''
		links_to_write = ''

		# Diploid organisms
		for step_id in range(0, int(param_dict['segsites']) * 2):
			#print('S\t'+ str(step_id + 1) + '\t' + ('0' if step_id < int(param_dict['segsites']) else '1'))
			steps_to_write += 'S\t'+ str(step_id + 1) + '\t' + ('0' if step_id < int(param_dict['segsites']) else '1') + '\n'

		#print('P\treference_haplotype\t' + ','.join([str(step_id) + '+' for step_id in range(0, int(param_dict['segsites']))]))
		paths_to_write = 'P\treference_haplotype\t' + ','.join([str(step_id) + '+' for step_id in range(1, int(param_dict['segsites']) + 1)]) + '\n'
			
		new_dict = {}
		for pos in range(0, int(param_dict['segsites'])):
			new_dict[pos] = [str(pos + 1), str(pos + 1 + 8)]

		#for pos, nodeIdREF_nodeIdALT_list in new_dict.items():
		#	print(pos, '-->', nodeIdREF_nodeIdALT_list)

		stuff_already_seen = set()
		for num_haplo, haplotype in enumerate(replicate_list[num_row:]):
			xxx_path_list = []
			for pos in range(0, len(haplotype)):
				node_in_pos = new_dict[pos][int(haplotype[pos])]
				xxx_path_list.append(node_in_pos + '+')

			#print('P\t' + 'individual_' + str(int(num_haplo/2) + 1) + '_haplotype_' + str(int(num_haplo%2) + 1) + '\t' + ','.join(xxx_path_list))
			paths_to_write += 'P\t' + 'individual_' + str(int(num_haplo/2) + 1) + '_haplotype_' + str(int(num_haplo%2) + 1) + '\t' + ','.join(xxx_path_list) + '\n'

			if haplotype not in stuff_already_seen: # if not already done, do it
				stuff_already_seen.add(haplotype)
				#print(haplotype)

				for pos in range(0, len(haplotype)-1):
					nodo_1 = new_dict[pos][int(haplotype[pos])]
					nodo_2 = new_dict[pos+1][int(haplotype[pos+1])]

					nodo_12 = nodo_1 + '_' + nodo_2
					if nodo_12 not in stuff_already_seen:
						stuff_already_seen.add(nodo_12) 

						#print('L\t' + str(nodo_1) + '\t+\t' + str(nodo_2) + '\t+\t0M')
						links_to_write += 'L\t' + str(nodo_1) + '\t+\t' + str(nodo_2) + '\t+\t0M\n'



		path_gfa = 'rep_' + str(num_rep) + '.gfa'
		with open(path_gfa, 'w') as fw:
			fw.write(steps_to_write + paths_to_write + links_to_write)
		
		print(path_gfa + ' written')


