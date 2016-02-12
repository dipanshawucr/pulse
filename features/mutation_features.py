# generate mutation features
from features_helpers import score_differences


def build_uniprot_to_index_to_disease(mutations_db):
    uniprot_to_index_to_disease = {}
    for line in mutations_db:
        try:
            tokens = line.split()
            uniprot = tokens[0]
            index = int(tokens[1])
            if uniprot in uniprot_to_index_to_disease:
                uniprot_to_index_to_disease[uniprot][index] = "*"
            else:
                uniprot_to_index_to_disease[uniprot] = {index: "*"}
        except ValueError:
            print "Cannot parse: " + line.strip()
    return uniprot_to_index_to_disease


def get_mutation_features(map_file, f_mutations_db_location, output_location):
    map_file_obj = open(map_file, 'r')
    mutations_db = open(f_mutations_db_location, 'r')
    write_to = open(output_location, 'w')
    uniprot_to_index_to_disease = build_uniprot_to_index_to_disease(mutations_db)
    for line in map_file_obj:
        tokens = line.split()
        asid = tokens[0].split("_")[0]
        prot = tokens[1]

        sstart = int(tokens[2])
        start = int(tokens[3])
        end = int(tokens[4])
        eend = int(tokens[5])

        if prot not in uniprot_to_index_to_disease:
            c1_count = 0
            a_count = 0
            c2_count = 0
            canonical_absolute = 0
        else:
            c1_count = score_differences(uniprot_to_index_to_disease, prot, sstart, start)
            a_count = score_differences(uniprot_to_index_to_disease, prot, start, end)
            c2_count = score_differences(uniprot_to_index_to_disease, prot, end, eend)
            prot_len = int(line.split("\t")[7].strip())
            canonical_absolute = score_differences(uniprot_to_index_to_disease, prot, 1, prot_len)

        print >> write_to, tokens[0] + "\t" + prot + "\t" + repr(c1_count) + "\t" + repr(a_count) + "\t" + repr(
            c2_count) + "\t" + repr(canonical_absolute)
    write_to.close()
