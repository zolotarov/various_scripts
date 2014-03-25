# pypy background2.py Athaliana_genome_cleanedup_utr.fasta 73638.46s user 11.29s system 99% cpu 20:28:53.91 total
# instead of 35 days
# pypy background_8mer_set.py Athaliana_genome_cleanedup_utr.fasta  17691.98s user 0.46s system 99% cpu 4:55:05.55 total

import itertools
from sys import argv
import operator

kmer_length = 8

all_kmers = [''.join(i) for i in itertools.product('ACGT', repeat = kmer_length)] #generate all kmers of the specified length

output_file = open('test8_set.bkgd', 'w')

def hamming_distance(kmer1, kmer2): # from http://code.activestate.com/recipes/499304-hamming-distance/
	return sum(itertools.imap(operator.ne, kmer1, kmer2))

def rev_comp(sequence):
	inseq =  "ACGT"
	outseq = "TGCA"
	sequence_unicode = u'%s' % sequence
	translation_table = dict(zip(map(ord, u'%s' % inseq), u'%s' % outseq))
	sequence_revcomp = sequence_unicode.translate(translation_table)[::-1]
	return str(sequence_revcomp)

def return_kmerred_sequences(file, reverese_complement = True):
	all_seq = list()
	with open(file, 'r') as input_file:
		sequence = ''
		for line in input_file:
			if ">" in line:
				all_seq.append(sequence)
				sequence = ''
				continue
			else:
				sequence += line.strip()
		all_seq.append(sequence)
	all_seq = all_seq[1:]
	joined = ''.join(all_seq)
	output_file.write('A\t%s\n' % joined.count('A'))
	output_file.write('C\t%s\n' % joined.count('C'))
	output_file.write('G\t%s\n' % joined.count('G'))
	output_file.write('T\t%s\n' % joined.count('T'))
	all_seq_reversed = [rev_comp(seq) for seq in all_seq]
	all_sequences_kmerred = list()
	if reverese_complement:
		for sequences in zip(all_seq, all_seq_reversed):
			kmers_in_sequence = [sequences[0][i:i + kmer_length] for i in range(len(sequences[0]) - kmer_length + 1)]
			kmers_in_sequence_rev = [sequences[1][i:i + kmer_length] for i in range(len(sequences[1]) - kmer_length + 1)]
			all_sequences_kmerred.append(set(kmers_in_sequence + kmers_in_sequence_rev))
	else:
		for sequence in all_seq:
			kmers_in_sequence = set([sequence[i:i + kmer_length] for i in range(len(sequence) - kmer_length + 1)])
			all_sequences_kmerred.append(kmers_in_sequence)
	return all_sequences_kmerred

def kmer_hd_dict(kmer):
	"""create a dictionary of lists of kmers and their corresponding 
		hamming distance to a given kmer

		kmer:	a kmer that will be compared to all other possible kmers
	"""
	hd_dict = {0:[], 1:[], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8:[]}
	for kmer2 in all_kmers:
		hd = hamming_distance(kmer, kmer2)
		hd_dict[hd].append(kmer2)
	for key in hd_dict:
		hd_dict[key] = set(hd_dict[key])
	return hd_dict

all_sequences_kmerred = return_kmerred_sequences(argv[1], True)

def counts(kmer):
	counts_dict = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0}
	kmer_hd = kmer_hd_dict(kmer)
	for list_of_kmers in all_sequences_kmerred:
		if kmer in list_of_kmers:
			counts_dict[0] += 1
			continue
		elif kmer_hd[1] & list_of_kmers:
			counts_dict[1] += 1
			continue
		elif kmer_hd[2] & list_of_kmers:
			counts_dict[2] += 1
			continue
		elif kmer_hd[3] & list_of_kmers:
			counts_dict[3] += 1
			continue
		elif kmer_hd[4] & list_of_kmers:
			counts_dict[4] += 1
			continue
		elif kmer_hd[5] & list_of_kmers:
			counts_dict[5] += 1
			continue
		elif kmer_hd[6] & list_of_kmers:
			counts_dict[6] += 1
			continue
		elif kmer_hd[7] & list_of_kmers:
			counts_dict[7] += 1
			continue
		elif kmer_hd[8] & list_of_kmers:
			counts_dict[8] += 1
			continue
	return counts_dict

for k in all_kmers:
	d = counts(k)
	e = ' '.join((str(d[i]) for i in d))
	output_file.write('%s\t%s\n' % (k, e))

output_file.close()
