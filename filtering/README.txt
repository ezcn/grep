
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ most_severe_consequence

Tipo di consequenza che è associata alla variante. 

Tipologie di consequenze che possono essere presenti nel file con il rispettivo ranking (che verrà riportato nella colonna soScore).
	{'transcript_ablation': 35, 'splice_acceptor_variant': 34, 'splice_donor_variant': 33, 'stop_gained': 32, 'frameshift_variant': 31, 'stop_lost': 30, 'start_lost': 29, 'transcript_amplification': 28, 'inframe_insertion': 27, 'inframe_deletion': 26, 'missense_variant': 25, 'protein_altering_variant': 24, 'splice_region_variant': 23, 'incomplete_terminal_codon_variant': 22, 'start_retained_variant': 21, 'stop_retained_variant': 20, 'synonymous_variant': 19, 'coding_sequence_variant': 18, 'mature_miRNA_variant': 17, '5_prime_UTR_variant': 16, '3_prime_UTR_variant': 15, 'non_coding_transcript_exon_variant': 14, 'intron_variant': 13, 'NMD_transcript_variant': 12, 'non_coding_transcript_variant': 11, 'upstream_gene_variant': 10, 'downstream_gene_variant': 9, 'TFBS_ablation': 8, 'TFBS_amplification': 7, 'TF_binding_site_variant': 6, 'regulatory_region_ablation': 5, 'regulatory_region_amplification': 4, 'feature_elongation': 3, 'regulatory_region_variant': 2, 'feature_truncation': 1, 'intergenic_variant': 0}

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ csqCount 

Conta dell'allele consequenza nell'individuo considerato.

Può assumere i valori : 0 (non presente), 1 (è presente in eterozigosi), 2 (è presente in omozigosi).

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ wCellCycle - wDDD - wEmbryoDev - wEssential - wLethal - wMisc

Viene controllato sei una variante genica appartiene o meno ad un lista dei geni fornita al programma. 

Può assumere i valori : 0 (non appartiene alla lista), 1 (appartiene alla lista).

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ rare 

Indica se la variante è rara o meno cercando nelle frequenze globali e valutando che in nessula delle popolazioni la frequenza superi una predeterminata threshold.

Può assumere i valori : 0 (supera la threshold in almeno una popolazione; non è rara), 0.5 (non possiedo i dati delle frequenze in nessuna popolazione; potrebbe essere rara o meno), 1 (non supera la threshold in nessuna delle popolazioni per le quali possiedo le frequenze; è rara).

 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ soScore

 E' un rank associato alla classificazione di VEP (variant effect predictor) in base alla gravità della consequenza.

 Assume il valore corrispettivo alla consequenza associata (Si può trovare del dizionario nel campo dedicato alla 'most_severe_consequence'

