def pbegend(pstart, pend, plength):
    phend = (pend + threshold / 3 > plength)
    phbeg = (pstart - threshold / 3 < 1)

    return [phbeg, phend]


def mbegend(nstart, nend, nlength, backwards):
    if backwards:
        temp = nstart
        nstart = nend
        nend = temp
    mhend = (nend + threshold > nlength)
    mhbeg = (nstart - threshold < 1)

    return [mhbeg, mhend]


# CONSTANTS
# TODO: Move these to preprocess_settings.json
threshold = 9  # nucleotides
pidentThreshold = 95.0
ignoreExclusion = False
dealWithNucleotide = False
lineNumber = 0


def filter_map1(blast_file, uniprot_fasta, isoform_fasta, rna_seq_output):
    read_from = open(blast_file, 'r')

    for line in read_from:

        (query_id,
         perc_identity,
         aln_length,
         query_start,
         query_end,
         subject_start,
         subject_end) = line.split('\t')

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

        if dealWithNucleotide:
            p_start /= 3.0
            p_end /= 3.0
            p_length /= 3.0

        nstart = int(query_start)
        nend = int(query_end)
        nlength = sum(nucleotide_lengths)

        pident = float(perc_identity)

        if query_id[0] == "E":
            if ignoreExclusion:
                continue
            nlength -= nucleotide_lengths[1]

        to_print = True

        backwards = (nstart > nend)

        [hit_p_beg, hit_p_end] = pbegend(p_start, p_end, p_length)
        [hit_n_beg, hit_n_end] = mbegend(nstart, nend, nlength, backwards)

        exist_beg = True
        exist_end = True

        if not hit_n_beg and not hit_p_beg:
            to_print = False
            exist_beg = False

        if not hit_n_end and not hit_p_end:
            to_print = False
            exist_end = False

        if not (exist_beg or exist_end):
            print line.strip() + " no_both"
        elif not exist_beg:
            if backwards:
                print line.strip() + " no_c2"
            else:
                print line.strip() + " no_c1"
        elif not exist_end:
            if backwards:
                print line.strip() + " no_c1"
            else:
                print line.strip() + " no_c2"
        else:
            hit_length = nend - nstart + 1
            if hit_length < 0:
                hit_length *= -1
                hit_length += 2

            if ((hit_length + threshold) / 3 < p_end - p_start + 1) or ((hit_length - threshold) / 3 > p_end - p_start + 1):
                to_print = False
                print line.strip() + " missing_chunks"

        if pident < pidentThreshold:
            to_print = False
            print line.strip() + " pident_failure"

        if query_id[0] == "E":
            # make sure A site is within
            if not backwards:
                if not (nstart <= nucleotide_lengths[0] <= nend):
                    to_print = False
            else:
                if not (nend <= nucleotide_lengths[0] <= nstart):
                    to_print = False

        if to_print:
            print >> rna_seq_output, line.strip()[0:len(line)]


def filter_map2(uniprot_fasta, isoform_fasta, filtermap_not_len_collapsed_output,
                                 filtermap_len_collapsed_output):

    # filters multiple uniprot hits by taking the best one in terms of evalue

    readFrom = open('./temp/rnaseq_huniprot_corrected_len.txt', 'r')
    writeTo = open('./temp/rnaseq_huniprot_corrected_len_collapsed.txt', 'w')

    oldID = ""

    for line in readFrom:
        tokens = line.split('\t')
        # PARSING ID
        asid = tokens[0]
        prot = tokens[1]
        # PARSING ID
        # id = asid + "//" + prot
        id = asid

        if oldID != id:
            print >> writeTo, line.strip()

        oldID = id

    ###
    # LOAD UNIPROT
    uniprot_ddbb = {}
    seq = ''
    uid = '|'
    SEP = '|'
    # print uniprot_fasta
    with open(uniprot_fasta, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                uniprot_ddbb[uid] = len(seq)
                # PARSING ID
                uid = line.split(SEP)[1]
                seq = ''
            # PARSING ID
            else:
                seq += line.strip()

    uniprot_ddbb[uid] = len(seq)
    ###
    # LOAD ISOFORM IN FASTA
    isoforms_ddbb = {}
    seq = ''
    uid = ''

    SEP = '>'
    with open(isoform_fasta, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                isoforms_ddbb[uid] = len(seq)
                # CHANGES THE LINE BELOW if is need it
                # PARSING ID
                uid = line.split(SEP)[1].strip()
                seq = ''
            # PARSING ID
            else:
                seq += line.strip()

    isoforms_ddbb[uid] = len(seq)

    # tags on spliced exon indices on query hits using pstart pend nstart and nend along with A exon positions
    # this information is obtained from the splice sequence to uniprot mapping file

    ##########################################################################
    ##########################################################################

    readFrom = open('./temp/rnaseq_huniprot_corrected_len_collapsed.txt', 'r')
    writeTo = open('uniprot_exon_indices_map.out', 'w')

    hasMapping = {}

    for line in readFrom:
        tokens = line.split('\t')
        # PARSING ID
        id = tokens[0]  # asid
        uniprot = tokens[1].split("|")[1]
        # PARSING ID
        if hasMapping.has_key(id):
            hasMapping[id].append(uniprot)
        else:
            hasMapping[id] = [uniprot]

    readFrom.close()

    ##########################################################################
    ##########################################################################

    readFrom = open('./temp/rnaseq_huniprot_corrected_len_collapsed.txt', 'r')

    for line in readFrom:
        tokens = line.split('\t')

        # PARSING ID
        id = tokens[0]  # asid
        uniprot = tokens[1].split("|")[1]
        # PARSING ID

        snucleotideLengths = tokens[0].split("-")[-1].split("=")
        nucleotideLengths = [int(snucleotideLengths[0]), int(snucleotideLengths[1]), int(snucleotideLengths[2])]

        pstart = int(tokens[8])
        pend = int(tokens[9])
        nstart = int(tokens[6])
        nend = int(tokens[7])

        if nstart < nend:  # these could be negative
            c1HitLength = nucleotideLengths[0] - nstart + 1
            c2HitLength = nucleotideLengths[2] - (sum(nucleotideLengths) - nend + 1)
        else:
            c1HitLength = nucleotideLengths[0] - nend + 1
            c2HitLength = nucleotideLengths[2] - (sum(nucleotideLengths) - nstart + 1)

        if c1HitLength < 0:
            c1HitLength = 0;
        if c2HitLength < 0:
            c2HitLength = 0;

        aStart = pstart + c1HitLength / 3
        aEnd = pend - c2HitLength / 3
        # PARSING ID
        # ID HANDELING change it if you need it
        transcript_isoform = id[2:].split('-')[0]
        # PARSING ID

        try:

            if aStart > 0 and aEnd > 0 and aEnd >= aStart and "-" not in uniprot:
                if id[0] == "I":
                    listProt = []
                    if hasMapping.has_key("E" + id[1:len(id)]):
                        listProt = hasMapping["E" + id[1:len(id)]]
                    print >> writeTo, id + "\t" + uniprot + "\t" + repr(pstart) + "\t" + repr(aStart) + "\t" + repr(
                        aEnd) + "\t" + repr(pend) + "\t" + repr(listProt) + "\t" + str(uniprot_ddbb[uniprot]) + "\t" + str(
                        isoforms_ddbb[transcript_isoform])
                else:
                    listProt = []
                    if hasMapping.has_key("I" + id[1:len(id)]):
                        listProt = hasMapping["I" + id[1:len(id)]]
                    print >> writeTo, id + "\t" + uniprot + "\t" + repr(pstart) + "\t" + repr(aStart) + "\t" + repr(
                        aStart) + "\t" + repr(pend) + "\t" + repr(listProt) + "\t" + str(
                        uniprot_ddbb[uniprot]) + "\t" + str(isoforms_ddbb[transcript_isoform])
            else:
                pass

        except KeyError, e:
            print '%s NOT FOUND ' % str(e)
            print 'Please check your Uniprot file or/and your poll of isoforms'
        # print id
