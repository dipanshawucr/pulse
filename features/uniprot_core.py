# reads uniprot core file and generates core features
from features_helpers import score_differences


def build_uniprot_to_index_to_core(sable_db_obj):
    uniprot_to_index_to_core = {}
    for line in sable_db_obj:
        tokens = line.split()
        try:
            # PARSING ID
            prot = tokens[0]
            index = int(tokens[1])
            core = tokens[2]
            # PARSING ID
            if uniprot_to_index_to_core.has_key(prot):
                uniprot_to_index_to_core[prot][index] = core
            else:
                uniprot_to_index_to_core[prot] = {index: core}
        except ValueError:
            print "Cannot parse: " + line[0:len(line) - 1]
    return uniprot_to_index_to_core


def get_sable_scores(map_file, f_sable_db_location, uniprot_core_output_location):
    map_file_obj = open(map_file, 'r')
    sable_db_obj = open(f_sable_db_location, 'r')
    write_to = open(uniprot_core_output_location, 'w')

    uniprot_to_index_to_core = build_uniprot_to_index_to_core(sable_db_obj)

    for line in map_file_obj:
        tokens = line.split()

        asid = tokens[0].split("_")[0]
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

        if prot in uniprot_to_index_to_core:
            c1_count = score_differences(uniprot_to_index_to_core, prot, sstart, start)
            a_count = score_differences(uniprot_to_index_to_core, prot, start, end)
            c2_count = score_differences(uniprot_to_index_to_core, prot, end, eend)
            prot_len = int(line.split("\t")[7].strip())
            canonical_absolute = score_differences(uniprot_to_index_to_core, prot, 1, prot_len)

        print >> write_to, tokens[0] + "\t" + prot + "\t" + repr(c1_count) + "\t" + repr(a_count) + "\t" + repr(
            c2_count) + "\t" + repr(canonical_absolute)
    write_to.close()