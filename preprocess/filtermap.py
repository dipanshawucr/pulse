def p_beg_end(p_start, pend, p_length):
    phend = (pend + THRESHOLD / 3 > p_length)
    phbeg = (p_start - THRESHOLD / 3 < 1)

    return [phbeg, phend]


def m_beg_end(n_start, n_end, n_length, backwards):
    if backwards:
        temp = n_start
        n_start = n_end
        n_end = temp
    mhend = (n_end + THRESHOLD > n_length)
    mhbeg = (n_start - THRESHOLD < 1)

    return [mhbeg, mhend]


# CONSTANTS
# TODO: Move these to preprocess_settings.json
THRESHOLD = 9  # nucleotides
PINDENT_THRESHOLD = 95.0
IGNORE_EXCLUSION = False
DEAL_WITH_NUCLEOTIDE = False


def filter_map1(blast_file, rna_seq_output):
    read_from = open(blast_file, 'r')
    write_to = open(rna_seq_output, 'w')
    for line in read_from:
        (query_id, subject_id, perc_identity, aln_length, mismatch_count, gap_open_count,
         query_start, query_end, subject_start, subject_end, e_val, bit_score) = line.split('\t')
        try:
            # PARSING ID
            snucleotide_lengths = query_id.split("-")[-1].split("=")
            nucleotide_lengths = [int(snucleotide_lengths[0]),
                                  int(snucleotide_lengths[1]),
                                  int(snucleotide_lengths[2])]
        except:
            print 'WARNING WITH' + query_id

        p_start = int(subject_start)  # /3.0
        p_end = int(subject_end)  # /3.0
        p_length = int(aln_length)  # /3.0

        if DEAL_WITH_NUCLEOTIDE:
            p_start /= 3.0
            p_end /= 3.0
            p_length /= 3.0

        n_start = int(query_start)
        n_end = int(query_end)
        n_length = sum(nucleotide_lengths)

        pident = float(perc_identity)

        if query_id[0] == "E":
            if IGNORE_EXCLUSION:
                continue
            n_length -= nucleotide_lengths[1]

        to_print = True

        backwards = (n_start > n_end)

        [hit_p_beg, hit_p_end] = p_beg_end(p_start, p_end, p_length)
        [hit_n_beg, hit_n_end] = m_beg_end(n_start, n_end, n_length, backwards)

        exist_beg = True
        exist_end = True

        if not hit_n_beg and not hit_p_beg:
            to_print = False
            exist_beg = False

        if not hit_n_end and not hit_p_end:
            to_print = False
            exist_end = False

        # if not (exist_beg or exist_end):
        #     print line.strip() + " no_both"
        # elif not exist_beg:
        #     if backwards:
        #         print line.strip() + " no_c2"
        #     else:
        #         print line.strip() + " no_c1"
        # elif not exist_end:
        #     if backwards:
        #         print line.strip() + " no_c1"
        #     else:
        #         print line.strip() + " no_c2"
        else:
            hit_length = n_end - n_start + 1
            if hit_length < 0:
                hit_length *= -1
                hit_length += 2

            if ((hit_length + THRESHOLD) / 3 < p_end - p_start + 1) or (
                    (hit_length - THRESHOLD) / 3 > p_end - p_start + 1):
                to_print = False
                # print line.strip() + " missing_chunks"

        if pident < PINDENT_THRESHOLD:
            to_print = False
            # print line.strip() + " pident_failure"

        if query_id[0] == "E":
            # make sure A site is within
            if not backwards:
                if not (n_start <= nucleotide_lengths[0] <= n_end):
                    to_print = False
            else:
                if not (n_end <= nucleotide_lengths[0] <= n_start):
                    to_print = False

        if to_print:
            print >> write_to, line.strip()[0:len(line)]


##########################################################################


def filter_map2(read_from, write_to):
    # filters multiple uniprot hits by taking the best one in terms of evalue
    read_from = open(read_from, 'r')
    write_to = open(write_to, 'w')

    old_id = ""

    for line in read_from:
        tokens = line.split('\t')
        # PARSING ID
        asid = tokens[0]
        _id = asid

        if old_id != _id:
            print >> write_to, line.strip()

        old_id = _id

##########################################################################


def filter_map3(uniprot_fasta):
    """
    Returns uniprot_ddbb.

    :param uniprot_fasta:
    :return:
    """
    uniprot_ddbb = {}
    seq = ''
    uid = '|'
    sep = '|'
    with open(uniprot_fasta, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                uniprot_ddbb[uid] = len(seq)
                # PARSING ID
                uid = line.split(sep)[1]
                seq = ''
            # PARSING ID
            else:
                seq += line.strip()

    uniprot_ddbb[uid] = len(seq)
    return uniprot_ddbb

##########################################################################


def filter_map4(isoform_fasta):
    """
    Load isoform in fasta.

    :param isoform_fasta:
    :return:
    """
    isoforms_ddbb = {}
    seq = ''
    uid = ''

    SEP = '>'
    with open(isoform_fasta, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                isoforms_ddbb[uid] = len(seq)
                # PARSING ID
                uid = line.split(SEP)[1].strip()
                seq = ''
            # PARSING ID
            else:
                seq += line.strip()

    isoforms_ddbb[uid] = len(seq)
    return isoforms_ddbb

##########################################################################


def filter_map5(read_from):
    """

    :param read_from:
    :return:
    """
    # tags on spliced exon indices on query hits using p_start pend n_start and n_end along with A exon positions
    # this information is obtained from the splice sequence to uniprot mapping file
    read_from = open(read_from, 'r')

    has_mapping = {}

    for line in read_from:
        tokens = line.split('\t')
        # PARSING ID
        id = tokens[0]  # asid
        uniprot = tokens[1].split("|")[1]
        # PARSING ID
        if has_mapping.has_key(id):
            has_mapping[id].append(uniprot)
        else:
            has_mapping[id] = [uniprot]

    read_from.close()
    return has_mapping

##########################################################################


def filter_map6(read_from, write_to, has_mapping, uniprot_ddbb, isoforms_ddbb):
    read_from = open(read_from, 'r')
    write_to = open(write_to, 'w')
    for line in read_from:
        tokens = line.split('\t')

        # PARSING ID
        _id = tokens[0]  # asid
        uniprot = tokens[1].split("|")[1]
        # PARSING ID

        snucleotide_lengths = tokens[0].split("-")[-1].split("=")
        nucleotide_lengths = [int(snucleotide_lengths[0]), int(snucleotide_lengths[1]), int(snucleotide_lengths[2])]

        p_start = int(tokens[8])
        p_end = int(tokens[9])
        n_start = int(tokens[6])
        n_end = int(tokens[7])

        if n_start < n_end:  # these could be negative
            c1_hit_length = nucleotide_lengths[0] - n_start + 1
            c2_hit_length = nucleotide_lengths[2] - (sum(nucleotide_lengths) - n_end + 1)
        else:
            c1_hit_length = nucleotide_lengths[0] - n_end + 1
            c2_hit_length = nucleotide_lengths[2] - (sum(nucleotide_lengths) - n_start + 1)

        if c1_hit_length < 0:
            c1_hit_length = 0
        if c2_hit_length < 0:
            c2_hit_length = 0

        a_start = p_start + c1_hit_length / 3
        a_end = p_end - c2_hit_length / 3

        # PARSING ID
        # ID HANDELING change it if you need it
        transcript_isoform = _id[2:].split('-')[0]
        # PARSING ID

        try:
            if (a_start and a_end) > 0 and a_end >= a_start and "-" not in uniprot:
                if _id[0] == "I":
                    list_prot = []
                    if ("E" + _id[1:len(_id)]) in has_mapping:
                        list_prot = has_mapping["E" + _id[1:len(_id)]]
                    print >> write_to, _id + "\t" + uniprot + "\t" + repr(p_start) + "\t" + repr(a_start) + \
                                       "\t" + repr(a_end) + "\t" + repr(p_end) + "\t" + repr(list_prot) + \
                                       "\t" + str(uniprot_ddbb[uniprot]) + "\t" + str(isoforms_ddbb[transcript_isoform])
                else:
                    list_prot = []
                    if ("I" + _id[1:len(_id)]) in has_mapping:
                        list_prot = has_mapping["I" + _id[1:len(_id)]]
                    print >> write_to, _id + "\t" + uniprot + "\t" + repr(p_start) + "\t" + repr(a_start) + \
                                       "\t" + repr(a_start) + "\t" + repr(p_end) + "\t" + repr(list_prot) + \
                                       "\t" + str(uniprot_ddbb[uniprot]) + "\t" + \
                                       str(isoforms_ddbb[transcript_isoform])
            else:
                pass

        except KeyError, e:
            print '%s NOT FOUND ' % str(e)
            print 'Please check your Uniprot file or/and your poll of isoforms'
            # print id
