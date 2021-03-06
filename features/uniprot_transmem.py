from features_helpers import score_differences


def get_transmembrane_region_features(uniprot_exon_indices_location, uniprot_tm_indices_db_location,
                                      output_file_location):
    """
    Reads uniprot TM indices DB and generates transmembrane scores.

    :param uniprot_exon_indices_location:
    :param uniprot_tm_indices_db_location:
    :param output_file_location:
    :return:
    """

    read_from = open(uniprot_exon_indices_location, 'r')
    write_to = open(output_file_location, 'w')
    uniprot_to_index_to_whatever = {}

    uniprot_tm_indices_db = open(uniprot_tm_indices_db_location, 'r')
    for line in uniprot_tm_indices_db:
        tokens = line.split("\t")
        try:
            uniprot = tokens[0].strip()
            start = int(tokens[2].strip())
            end = int(tokens[3].strip())
            for i in range(start, end + 1):
                if uniprot_to_index_to_whatever.has_key(uniprot):
                    uniprot_to_index_to_whatever[uniprot][i] = "*"
                else:
                    uniprot_to_index_to_whatever[uniprot] = {i: "*"}
        except ValueError:
            print "Cannot parse: " + line[0:len(line) - 1]

    for line in read_from:
        tokens = line.split()
        # PARSING ID
        asid = tokens[0].split("_")[0]
        prot = tokens[1]

        # PARSING ID
        sstart = int(tokens[2])
        start = int(tokens[3])
        end = int(tokens[4])
        eend = int(tokens[5])

        c1_count = 0
        a_count = 0
        c2_count = 0
        canonical_absolute = 0

        if not uniprot_to_index_to_whatever.has_key(prot):
            c1_count = 0
            a_count = 0
            c2_count = 0
            canonical_absolute = 0
            other_absolute = 0
        else:
            c1_count = score_differences(uniprot_to_index_to_whatever, prot, sstart, start)
            a_count = score_differences(uniprot_to_index_to_whatever, prot, start, end)
            c2_count = score_differences(uniprot_to_index_to_whatever, prot, end, eend)
            prot_len = int(line.split("\t")[7].strip())
            canonical_absolute = score_differences(uniprot_to_index_to_whatever, prot, 1, prot_len)
        print >> write_to, tokens[0] + "\t" + prot + "\t" + repr(c1_count) + "\t" + repr(a_count) + "\t" + repr(
            c2_count) + "\t" + repr(canonical_absolute)
