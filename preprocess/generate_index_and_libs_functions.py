def translate_backwards(seq, query_start, query_end):
    tn_seq = seq.complement()
    start, end = reverse(seq, query_start, query_end)
    n_seq = str(tn_seq)[::-1]

    bases = ['T', 'C', 'A', 'G']
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids))
    p_seq = ""
    cdna_seq = []

    # the start point, is where the alignment begins, not the protein, that give us the read frame
    # get the seq from this point until I found a stop codon
    i = start
    while True:
        if i <= 0:
            break
        codon = n_seq[i - 3:i]
        i -= 3
        if len(codon) != 3:
            break
        amino_acid = codon_table[codon]
        if amino_acid == "*":
            break
        p_seq += amino_acid
        cdna_seq.append(codon)
    p_seq = p_seq[::-1]
    cdna_seq = cdna_seq[::-1]
    i = start
    while True:
        codon = n_seq[i:i + 3]
        i += 3
        if len(codon) != 3:
            break
        amino_acid = codon_table[codon]
        if amino_acid == "*":
            break
        p_seq += amino_acid
        cdna_seq.append(codon)

    # Remove all the AA of the begin of the seq that are not a M
    for a in range(len(p_seq)):
        if p_seq[a] != 'M':
            pass
        else:
            p_seq = p_seq[a:]
            cdna_seq = cdna_seq[a:]
            break
    cdna_seq = ''.join(cdna_seq)
    return p_seq, cdna_seq, start, end


def translate(n_seq, start):
    bases = ['T', 'C', 'A', 'G']
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids))
    p_seq = ""
    cdna_seq = []

    # the start point, is where the alignment begins, not the protein, that give us the read frame
    # get the seq from this point until I found a stop codon

    i = start - 1
    while True:
        codon = n_seq[i - 3:i]
        i -= 3
        if len(codon) != 3:
            break
        amino_acid = codon_table[codon]
        if amino_acid == "*":
            break
        p_seq += amino_acid
        cdna_seq.append(codon)
    p_seq = p_seq[::-1]
    cdna_seq = cdna_seq[::-1]

    i = start - 1
    while True:
        codon = n_seq[i:i + 3]
        i += 3
        if len(codon) != 3:
            break
        amino_acid = codon_table[codon]
        if amino_acid == "*":
            break
        p_seq += amino_acid
        cdna_seq.append(codon)

    # Remove all the AA of the begin of the seq that are not a M
    for a in range(len(p_seq)):
        if p_seq[a] == 'M':
            p_seq = p_seq[a:]
            cdna_seq = cdna_seq[a:]
            break
        else:
            pass
    cdna_seq = ''.join(cdna_seq)
    return p_seq, cdna_seq


def reverse(seq, start, end):
    r_start = len(seq) - (int(start) + 1)
    r_end = len(seq) - int(end)
    return r_start, r_end


def get_n_seq(query_id, complete_transcripts_location):
    nucleo_transcript_seqs = open(complete_transcripts_location, 'r')
    read = False
    qid = query_id.split('\"')[1]

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
