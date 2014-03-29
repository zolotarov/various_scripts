import itertools
import operator
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="input FASTA file", required=True)
parser.add_argument(
    "-o", "--output", help="output background distribution",
    required=True)
parser.add_argument(
    "-l", "--length", help="length of the kmer to analyze",
    default='6', choices=['6', '8'])
parser.add_argument(
    "-s", "--strand", help="strands to analyze, both or just forward",
    default="both", choices=['both', 'forward'])
args = vars(parser.parse_args())
strand = args["strand"]
kmer_length = int(args["length"])

# generate all kmers of the specified length
all_kmers = [''.join(i) for i in itertools.product('ACGT',
             repeat=kmer_length)]

output_file = open(args["output"], 'w')


def hamming_distance(kmer1, kmer2):
    # from http://code.activestate.com/recipes/499304-hamming-distance/
    return sum(itertools.imap(operator.ne, kmer1, kmer2))


def rev_comp(sequence):
    """return a reverse complement of a sequence
        sequence:   a FASTA sequence in capitals without ambiguous nucleotides
    """
    inseq = "ACGT"
    outseq = "TGCA"
    sequence_unicode = u'%s' % sequence
    translation_table = dict(zip(map(ord, u'%s' % inseq), u'%s' % outseq))
    sequence_revcomp = sequence_unicode.translate(translation_table)[::-1]
    return str(sequence_revcomp)


def return_kmerred_sequences(file, strand):
    """take all sequences and convert them into list of lists of kmers

    file:   the input filename
    reverse_complement: True or False, whether to find smalles HD
    in forward orientation or both
    """
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
    if strand == 'both':
        for sequences in zip(all_seq, all_seq_reversed):
            kmers_in_sequence = [
                sequences[0][i:i + kmer_length] for i in
                range(len(sequences[0]) - kmer_length + 1)]
            kmers_in_sequence_rev = [
                sequences[1][i:i + kmer_length] for i in
                range(len(sequences[1]) - kmer_length + 1)]
            all_sequences_kmerred.append(set(kmers_in_sequence
                                         + kmers_in_sequence_rev))
    else:
        for sequence in all_seq:
            kmers_in_sequence = set(
                [sequence[i:i + kmer_length] for i in
                    range(len(sequence) - kmer_length + 1)])
            all_sequences_kmerred.append(kmers_in_sequence)
    return all_sequences_kmerred


def kmer_hd_dict(kmer):
    """create a dictionary of lists of kmers and their corresponding
        hamming distance to a given kmer

        kmer:   a kmer that will be compared to all other possible kmers
    """
    hd_dict = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: []}
    for kmer2 in all_kmers:
        hd = hamming_distance(kmer, kmer2)
        hd_dict[hd].append(kmer2)
    for key in hd_dict:
        hd_dict[key] = set(hd_dict[key])
    return hd_dict

all_sequences_kmerred = return_kmerred_sequences(args["input"], strand)


def counts(kmer):
    """create a dictionary of all SMDs (substring minimal distance), i.e.
        sum of the number of sequences containing a minimum Hamming distance
        to a kmer

        kmer: a kmer for which SMDs are calculated for all sequences in the
              background set
    """
    if kmer_length == 6:
        counts_dict = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
    elif kmer_length == 8:
        counts_dict = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0}

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
