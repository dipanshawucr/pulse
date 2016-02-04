from features_helpers import score_differences


def get_uniprot_elm_features(uniprot_exon_indices_location, uniprot_elm_db_location, output_location):
    """
    Reads uniprot ELM file and generates ELM features.

    :param uniprot_exon_indices_location:
    :param uniprot_elm_db_location:
    :param output_location:
    :return:
    """

    read_from = open(uniprot_elm_db_location, 'r')
    uniprot_exon_indices = open(uniprot_exon_indices_location, 'r')
    write_to = open(output_location, 'w')
    uniprot_to_index_to_elm = {}

    for line in read_from:
        tokens = line.split("\t")
        try:
            uniprot = tokens[0].strip()
            start = int(tokens[2].strip())
            end = int(tokens[3].strip())
            for i in range(start, end + 1):
                if uniprot in uniprot_to_index_to_elm:
                    uniprot_to_index_to_elm[uniprot][i] = "*"
                else:
                    uniprot_to_index_to_elm[uniprot] = {i: "*"}
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

        if not uniprot_to_index_to_elm.has_key(prot):
            c1_count = 0
            a_count = 0
            c2_count = 0
            canonical_absolute = 0
            other_absolute = 0
        else:
            c1_count = score_differences(uniprot_to_index_to_elm, prot, sstart, start)
            a_count = score_differences(uniprot_to_index_to_elm, prot, start, end)
            c2_count = score_differences(uniprot_to_index_to_elm, prot, end, eend)
            protLen = int(line.split("\t")[7].strip())
            canonical_absolute = score_differences(uniprot_to_index_to_elm, prot, 1, protLen)
        print >> write_to, tokens[0] + "\t" + prot + "\t" + repr(c1_count) + "\t" + repr(a_count) + "\t" + repr(
            c2_count) + "\t" + repr(canonical_absolute)
