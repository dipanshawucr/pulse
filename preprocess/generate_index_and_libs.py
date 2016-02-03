from Bio.Seq import Seq
import sys


def translate_backwards(seq, queryStart, queryEnd):
    tn_seq = seq.complement()
    start, end = reverse(seq, queryStart, queryEnd)
    nseq = str(tn_seq)[::-1]

    # print nseq
    bases = ['T', 'C', 'A', 'G']
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids))
    pseq = ""
    cdna_seq = []

    # the start point, is where the aligment begins, not the protein, that give us the read frame
    # get the seq from this point until I found a stop codon
    i = start
    while True:
        if i <= 0:
            # print 'entra'
            break
        codon = nseq[i - 3:i]
        # print codon
        i -= 3
        # print len(codon),codon
        if len(codon) != 3:
            # print 'codon small'
            break

        aminoAcid = codon_table[codon]

        if aminoAcid == "*":
            # print 'amino acid *'
            break

        pseq += aminoAcid
        cdna_seq.append(codon)
        # print aminoAcid
        # print pseq

    # print len(pseq)
    # reverse
    pseq = pseq[::-1]
    cdna_seq = cdna_seq[::-1]
    # print len(pseq)
    # print 'next'

    i = start
    while True:
        # print 'entro en el segundo bucle'
        codon = nseq[i:i + 3]
        i += 3
        if len(codon) != 3:
            break
        aminoAcid = codon_table[codon]
        if aminoAcid == "*":
            break
        pseq += aminoAcid
        cdna_seq.append(codon)

    # Remove all the AA of the begin of the seq that are not a M
    for a in range(len(pseq)):
        # print pseq[a]
        if pseq[a] != 'M':
            pass
        else:
            pseq = pseq[a:]
            cdna_seq = cdna_seq[a:]
            break
    cdna_seq = ''.join(cdna_seq)
    return pseq, cdna_seq, start, end


def translate(nseq, start):
    bases = ['T', 'C', 'A', 'G']
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids))
    pseq = ""
    cdna_seq = []

    # the start point, is where the aligment begins, not the protein, that give us the read frame
    # get the seq from this point until I found a stop codon

    i = start - 1
    while True:
        codon = nseq[i - 3:i]
        i -= 3
        if len(codon) != 3:
            break
        aminoAcid = codon_table[codon]
        if aminoAcid == "*":
            break
        pseq += aminoAcid
        cdna_seq.append(codon)
    pseq = pseq[::-1]
    cdna_seq = cdna_seq[::-1]

    i = start - 1
    while True:
        codon = nseq[i:i + 3]
        i += 3
        if len(codon) != 3:
            break
        aminoAcid = codon_table[codon]
        if aminoAcid == "*":
            break
        pseq += aminoAcid
        cdna_seq.append(codon)

    # Remove all the AA of the begin of the seq that are not a M
    for a in range(len(pseq)):
        # print pseq[a]
        if pseq[a] == 'M':
            pseq = pseq[a:]
            cdna_seq = cdna_seq[a:]
            break
        else:
            pass
    cdna_seq = ''.join(cdna_seq)
    return pseq, cdna_seq


def reverse(seq, start, end):
    rstart = len(seq) - (int(start) + 1)
    rend = len(seq) - int(end)
    return rstart, rend


def get_nseq(query_id):
    # full transcript file
    nucleo_transcript_seqs = open('complete_transcripts.fasta', 'r')
    read = False
    qid = query_id.split('\"')[1]
    # qid = '>'+qid

    for line in nucleo_transcript_seqs:
        if line[0] == '>':
            ddbb_id = line.split('\"')[1]
            if ddbb_id == qid:
                read = True
        elif read:
            return line.strip()
        else:
            pass
    return


try:
    blastx_file = sys.argv[1]
    blastx_output = open(blastx_file, 'r')
except:
    print 'Fatal Error, blast index must be provided'
    print 'Remember to provide a file filtered, like rnaseq_huniprot_corrected_len_collapsed.txt'
    sys.exit()
try:
    open('complete_transcripts.fasta', 'r')
except:
    print 'Fatal Error, Transcript sequence must be provided'
    print 'nseq_isoforms.fas, not found'
    sys.exit()

index = 'uniprot_exon_indices_map.out'
filter_index = {}
with open(index, 'r') as input_file:
    for line in input_file:
        data = line.split('\t')
        filter_index[data[0]] = 1

cdna_output = open('cdna_transcripts.fasta', 'w')
pseq_output = open('pSeq_isoforms.fas', 'w')
# blastx_output = open('../preprocess/rnaseq_huniprot_corrected_len_collapsed.txt', 'r')
# blastx_output = open('./Data_Ath_Pulse_run/test_transcripts_len.out', 'r')
# nucleo_transcript_seqs =open('./Data_Ath_Pulse_run/aradobsis_transcripts.fas','r')

for blastx_hit in blastx_output:

    try:
        (queryId, subjectId, percIdentity, alnLength, mismatchCount, gapOpenCount, queryStart, queryEnd, subjectStart,
         subjectEnd, eVal, bitScore) = blastx_hit.split('\t')
        if filter_index.has_key(queryId):

            n_seq = Seq(get_nseq(queryId))

            if int(queryStart) > int(queryEnd):

                pseq, cdna_seq, s, e = translate_backwards(n_seq, int(queryStart), int(queryEnd))
                # Save AA seq
                print >> pseq_output, '>' + str(queryId)
                print >> pseq_output, pseq
                # Save CDNA seq
                print >> cdna_output, '>' + str(queryId)
                print >> cdna_output, cdna_seq

            else:
                tn_seq = str(n_seq)
                pseq, cdna_seq = translate(tn_seq, int(queryStart))

                print >> pseq_output, '>' + str(queryId)
                print >> pseq_output, pseq
                # Save CDNA seq
                print >> cdna_output, '>' + str(queryId)
                print >> cdna_output, cdna_seq

        else:
            pass
    except Exception, e:
        print e

pseq_output.close()
cdna_output.close()


# load putative isoforms lenght
dict_isoforms_len = {}

for line in open('pSeq_isoforms.fas', 'r'):
    if line[0] == '>':

        # id = line.split('\"')[1]
        sid = line[1:].strip()

        dict_isoforms_len[sid] = 0
    else:
        dict_isoforms_len[sid] = len(line.strip())
# load canonical form lenght
# load ids
print len(dict_isoforms_len)
# print dict_isoforms_len['I."ENST00000542353";12_7898927_7899019-111=93=161']
dict_canonical_len = {}

for line in open('uniprot_exon_indices_map.out', 'r'):
    tab_data = line.split('\t')
    dict_canonical_len[tab_data[1]] = 0

seq_len = 0
write = False
for line in open('/home/wonjunetai/Documents/KimLab/WJ_09_292015_test/input/uniprot_sprot.fasta', 'r'):
    if line[0] == '>':
        if write:
            dict_canonical_len[protein_id] = seq_len
            seq_len = 0
            write = False

        protein_id = line.split('|')[1]

        if dict_canonical_len.has_key(protein_id):
            write = True
            seq_len = 0

    else:
        seq_len += len(line.strip())

# feacth length
output = open('uniprot_exon_indices.txt', 'w')
for line in open('uniprot_exon_indices.txt', 'r'):
    tab_data = line.split('\t')
    canonical_id = tab_data[1]
    transcript_id = tab_data[0].split('\"')[1]
    transcript_id = tab_data[0].strip()
    # print dict_canonical_len[canonical_id],dict_isoforms_len[transcript_id]
    print >> output, line.strip() + '\t' + str(dict_canonical_len[canonical_id]) + '\t' + str(
        dict_isoforms_len[transcript_id])
