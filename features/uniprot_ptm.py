from features_helpers import score_differences


def get_postranscriptional_modification_features(uniprot_exon_indices_location, uniprot_ptm_db_location,
                                                 output_file_location):
    """
    Reads uniprot PTM file and generates post translational modification site features.

    :param uniprot_exon_indices_location:
    :param uniprot_ptm_db_location:
    :param output_file_location:
    :return:
    """

    read_from = open(uniprot_ptm_db_location, 'r')
    uniprot_exon_indices = open(uniprot_exon_indices_location, 'r')
    write_to = open(output_file_location, 'w')
    uniprot_to_index_to_ptm = {}

    for line in read_from:
        tokens = line.split()
        try:
            uniprot = tokens[0]
            index = int(tokens[1])
            ptm = tokens[3]
            if uniprot_to_index_to_ptm.has_key(uniprot):
                uniprot_to_index_to_ptm[uniprot][index] = "*"
            else:
                uniprot_to_index_to_ptm[uniprot] = {index: "*"}
        except ValueError:
            print "Cannot parse: " + line[0:len(line) - 1]

    for line in uniprot_exon_indices:
        tokens = line.split()

        asid = tokens[0].split("_")[0]
        prot = tokens[1]

        sstart = int(tokens[2])
        start = int(tokens[3])
        end = int(tokens[4])
        eend = int(tokens[5])

        c1_count = 0
        a_count = 0
        c2_count = 0
        canonical_absolute = 0

        if not uniprot_to_index_to_ptm.has_key(prot):
            c1_count = 0
            a_count = 0
            c2_count = 0
            canonical_absolute = 0
            other_absolute = 0
        else:
            c1_count = score_differences(uniprot_to_index_to_ptm, prot, sstart, start)
            a_count = score_differences(uniprot_to_index_to_ptm, prot, start, end)
            c2_count = score_differences(uniprot_to_index_to_ptm, prot, end, eend)
            prot_len = int(line.split("\t")[7].strip())
            canonical_absolute = score_differences(uniprot_to_index_to_ptm, prot, 1, prot_len)
        print >> write_to, tokens[0] + "\t" + prot + "\t" + repr(c1_count) + "\t" + repr(a_count) + "\t" + repr(
            c2_count) + "\t" + repr(canonical_absolute)
