# reads uniprot disorder file and generates disorder features
# do absolute disorder as well
from features_helpers import score_differences, link_db


def fetch_feature(feature, anchors_map, path_ddbb):
    (c, db) = link_db(path_ddbb)
    dict_feature = {}
    for anchor in anchors_map:
        # print anchor
        query = c.execute(""" SELECT * FROM %s WHERE id='%s' """ % (feature, anchor.strip()))
        # query = c.execute("""SELECT * FROM uniprot_ptms WHERE uniprot_id = "Q8IZP0" """)
        # Row_id Uniprot_id postion Value ['*'] in pfam ['*',enzymatic bool]
        for row in query:
            uniprot_id = row[0]
            index = int(row[1])
            value = row[2]

            if dict_feature.has_key(uniprot_id):
                dict_feature[uniprot_id][index] = value
            else:
                dict_feature[uniprot_id] = {index: value}
    c.close()
    return dict_feature


def build_anchors_map(map_file, anchors_map=list()):
    uniprot_exon_indices = open(map_file, 'r')
    for line in uniprot_exon_indices:
        tokens = line.split('\t')
        prot = tokens[1]
        anchors_map.append(prot)
    uniprot_exon_indices.close()
    return anchors_map


def build_uniprot_to_index_to_disorder(iupred_output, index_to_disorder):
    iupred_output_file = open(iupred_output, 'r')
    for line in iupred_output_file:
        tokens = line.split('\t')
        if len(tokens) > 1:
            try:
                prot = tokens[0]
                index = int(tokens[1])
                disordered = tokens[2].strip()
                if prot in index_to_disorder:
                    index_to_disorder[prot][index] = disordered
                else:
                    index_to_disorder[prot] = {index: disordered}
            except ValueError:
                print "Cannot parse: " + line[0:len(line) - 1]
    iupred_output_file.close()
    return index_to_disorder


def get_uniprot_disorder_features(pulse_path, map_file, iupred_output,
                                  canonical_db_location, disorder_read_out_location):

    anchors_map = build_anchors_map(map_file)
    uniprot_to_index_to_disorder = fetch_feature('uniprot_iupred', anchors_map, canonical_db_location)
    uniprot_to_index_to_disorder = build_uniprot_to_index_to_disorder(iupred_output, uniprot_to_index_to_disorder)

    write_to = open(disorder_read_out_location, 'w')
    uniprot_exon_indices = open(map_file, 'r')
    for line in uniprot_exon_indices:
        tokens = line.split('\t')
        asid = tokens[0]  # .split("_")[0]
        prot = tokens[1]
        sstart = int(tokens[2])
        start = int(tokens[3])
        end = int(tokens[4])
        eend = int(tokens[5])
        rough_a_length = int(int(tokens[0].split("_")[-1].split("=")[1]) / 3)
        if asid[0] == "I":
            rough_a_length = 0
        c1_count = 0
        a_count = 0
        c2_count = 0
        canonical_absolute = 0
        other_absolute = 0
        other_a = 0
        if prot not in uniprot_to_index_to_disorder:
            print "Protein not in uniprot_to_index_to_disorder: ", prot
            c1_count = 0
            a_count = 0
            c2_count = 0
            canonical_absolute = 0
            other_absolute = 0
            other_a = 0
        else:
            c1_count = score_differences(uniprot_to_index_to_disorder, prot, sstart, start)
            a_count = score_differences(uniprot_to_index_to_disorder, prot, start, end)
            c2_count = score_differences(uniprot_to_index_to_disorder, prot, end, eend)
            prot_len = int(line.split("\t")[7].strip())
            other_prot = asid  # or token[0] Customize this line if is need it for test [asid+'//'+prot+'-ST']
            other_prot_len = int(line.split("\t")[8].strip())
            canonical_absolute = score_differences(uniprot_to_index_to_disorder, prot, 1, prot_len)
            other_absolute = score_differences(uniprot_to_index_to_disorder, other_prot, 1, other_prot_len)
            other_a_end = start + rough_a_length

            if other_a_end > other_prot_len:
                other_a_end = other_prot_len
            otherA = score_differences(uniprot_to_index_to_disorder, other_prot, start, other_a_end)
        print >> write_to, tokens[0] + "\t" + prot + "\t" + repr(c1_count) + "\t" + repr(a_count) + "\t" + repr(
            c2_count) + "\t" + repr(canonical_absolute) + "\t" + repr(other_absolute) + "\t" + repr(other_a)
    write_to.close()
