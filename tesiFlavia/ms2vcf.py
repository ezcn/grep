import os

replicate_to_pos_locus_to_sample_to_genotype_dict = {}

with open('/home/flavia/Desktop/ms_to_vcf/simtwopop.ms') as f:
	# Skit the first 2 rows
	for _ in range(2):
		next(f)

	# Fast implementation (to improve for huge simulations)
	for num_rep, replicate in enumerate(f.read().split('//')[1:]):
		replicate_to_pos_locus_to_sample_to_genotype_dict[num_rep] = {}

		replicate_list = replicate.strip().split('\n')
		segsites_num = int(replicate_list[0].split('segsites: ')[1])
		positions_list = replicate_list[1].split('positions: ')[1].strip().split(' ')

		#print('segsites_num:', segsites_num)
		#print('positions_list:', positions_list)

		num_sample = 0
		haplo_iterator = iter(replicate_list[2:])
		for haplotype1 in haplo_iterator:
			haplotype2 = next(haplo_iterator)
			#print(haplotype1, haplotype2)

			for pos, locus_haplo1, locus_haplo2 in zip(positions_list, haplotype1, haplotype2):
				tmp = int(locus_haplo1) + int(locus_haplo2)
				if tmp == 0:                        #if is one 0 p
					genotype = '0/0'
				elif tmp == 1:                      #if 1
					genotype = '0/1'
				else:
					genotype = '1/1'

				if pos not in replicate_to_pos_locus_to_sample_to_genotype_dict[num_rep].keys():
					replicate_to_pos_locus_to_sample_to_genotype_dict[num_rep][pos] = {}
				
				replicate_to_pos_locus_to_sample_to_genotype_dict[num_rep][pos][num_sample] = genotype
				#print(pos, locus_haplo1, locus_haplo2, genotype, 'sample' + str(num_sample))

			num_sample += 1

from datetime import date


for num_rep, pos_locus_to_sample_to_genotype_dict in replicate_to_pos_locus_to_sample_to_genotype_dict.items():
	with open('ms_rep{}.vcf'.format(num_rep + 1), 'w') as fw:
		fw.write('##fileformat=VCFv4.2\n')
		fw.write('##fileDate={}\n'.format(date.today().strftime("%Y%m%d")))
		fw.write('##reference=unknown.fa\n')
		fw.write('##commandline=ms2vcf.py --input-ms xxx.ms\n')
		fw.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
		fw.write('\t'.join(
			['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + ['Sample_{}'.format(i + 1) for i in range(num_sample)]
		) + '\n')
		for pos, sample_to_genotype_dict in pos_locus_to_sample_to_genotype_dict.items():
			vcf_row_list = ['Unkwn_rep' + str(num_rep + 1), pos, '.', '.', '.', '.', '.', '.', 'GT']
			for n_sample in range(num_sample):
				vcf_row_list += [sample_to_genotype_dict[n_sample]]

			fw.write('\t'.join(vcf_row_list) + '\n')
