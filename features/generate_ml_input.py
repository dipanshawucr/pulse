# Combines all our feature tables into one huge feature table ready to be read by R
# Also generates some simple features not found in other feature tables


def read_features(map, stream, duality):
    for line in stream:
        tokens = line.split()
        if duality:
            key = tokens[0] + "//" + tokens[1]
            rest = tokens[2:len(tokens)]
            map[key] = rest
        else:
            key = tokens[0]
            rest = tokens[1:len(tokens)]
            map[key] = rest


def obtain_feature(mapp, key):
    try:
        length = len(mapp[mapp.keys()[0]])
    except TypeError:
        length = 1
    new_list = ["0"] * length
    if key in mapp:
        new_list = mapp[key]
    return new_list


def convert_list(input_list):
    try:
        string_r = ""
        for i in input_list:
            string_r += i + ", "
        return string_r[0:len(string_r) - 2]
    except TypeError:
        return repr(input_list)


def generate_machine_learning_matrix(map_file, core_read, degree_read, disorder_read, domain_read, elm_read,
                                     event_con_read, mutation_read, ptm_read, sequence_con_read, transmem_read,
                                     write_all_output, write_all_names_output):
    map_file_obj = open(map_file, 'r')

    read_core = open(core_read, 'r')
    read_degree = open(degree_read, 'r')
    read_disorder = open(disorder_read, 'r')
    read_domain = open(domain_read, 'r')
    read_elm = open(elm_read, 'r')
    read_event_con = open(event_con_read, 'r')
    read_mutation = open(mutation_read, 'r')
    read_ptm = open(ptm_read, 'r')
    read_sequence_con = open(sequence_con_read, 'r')
    read_transmem = open(transmem_read, 'r')

    # This feature must be make ad hoc for our data (Event classification)
    # Please, read the supplementary Material
    # or just ignore this feature, notice the global score distribution will be affected
    # readTS = open(PulsePATH+'/src_features/param_files/ML_AS_Event_Classification.inp', 'r')

    write_all = open(write_all_output, 'w')
    write_all_names = open(write_all_names_output, 'w')

    # dual key features
    core_read = {}
    disorder_read = {}
    domain_read = {}
    elm_read = {}
    mutation_read = {}
    ptm_read = {}
    transmem_read = {}

    # single key features
    degree_read = {}
    sequenceCon_read = {}
    eventCon_read = {}
    ts_read = {}

    read_features(core_read, read_core, True)
    read_features(disorder_read, read_disorder, True)
    read_features(domain_read, read_domain, True)
    read_features(elm_read, read_elm, True)
    read_features(mutation_read, read_mutation, True)
    read_features(ptm_read, read_ptm, True)
    read_features(transmem_read, read_transmem, True)
    read_features(degree_read, read_degree, False)
    read_features(sequenceCon_read, read_sequence_con, False)
    read_features(eventCon_read, read_event_con, False)

    # for line in readTS:
    #     tokens = line.split()
    #     asid = tokens[0]
    #     ts = 0
    #     if tokens[2] == "TS":
    #         ts = 1
    #     ts_read[asid] = ts

    output_positiveI = []
    output_uniprotI = []
    output_restI = []
    output_positiveE = []
    output_uniprotE = []
    output_restE = []

    for line in map_file_obj:
        tokens = line.split()
        asid = tokens[0]
        prot = tokens[1]
        sstart = int(tokens[2])
        start = int(tokens[3])
        end = int(tokens[4])
        eend = int(tokens[5])
        canonical_length = int(line.split("\t")[7].strip())
        other_length = int(line.split("\t")[8].strip())

        temp_line = line.split("\t")[6].strip().replace("[", "").replace("]", "").replace("'", "").replace(" ", "")

        list_other = []
        if temp_line != "":
            list_other = temp_line.split(",")

        # check the list to see if it's uniprot (not empty), or positive and add the
        # line to print which we get eventually.
        # at the end print these in order.

        dual_key = asid + "//" + prot
        single_key = asid

        # obtain feature line from various feature files.
        core_features = obtain_feature(core_read, dual_key)
        disorder_features = obtain_feature(disorder_read, dual_key)
        domain_features = obtain_feature(domain_read, dual_key)
        elm_features = obtain_feature(elm_read, dual_key)
        mutation_features = obtain_feature(mutation_read, dual_key)
        ptm_features = obtain_feature(ptm_read, dual_key)
        transmem_features = obtain_feature(transmem_read, dual_key)

        degree_features = obtain_feature(degree_read, single_key)
        sequence_con_features = obtain_feature(sequenceCon_read, single_key)
        event_con_features = obtain_feature(eventCon_read, single_key)
        ts_features = obtain_feature(ts_read, single_key)

        len_diff = canonical_length - other_length
        n_len_diff = (1.0 * len_diff) / canonical_length
        len_after = other_length - start
        frameshift = 1

        # This is ready for the default format label splicing plus lenght of the exons -> X.LABEL-###_33=33=33
        if int(tokens[0].split("_")[-1].split("=")[1]) % 3 == 0:
            frameshift = 0

        length_a_exon = int(tokens[0].split("-")[-1].split("=")[1])

        feature_line = repr(frameshift) + ", " + repr(length_a_exon) + ", " + repr(len_diff) + ", " + \
                       repr(n_len_diff) + ", " + repr(len_after) + ", " + convert_list(core_features) + ", " + \
                       convert_list(disorder_features) + ", " + convert_list(domain_features) + ", " + \
                       convert_list(elm_features) + ", " + convert_list(mutation_features) + ", " + \
                       convert_list(ptm_features) + ", " + convert_list(transmem_features) + ", " + \
                       convert_list(degree_features) + ", " + convert_list(sequence_con_features) + ", " + \
                       convert_list(event_con_features) + ", " + convert_list(ts_features)

        if asid[0] == "I":
            output_positiveI.append([dual_key, feature_line])

        if asid[0] == "E":
            output_positiveE.append([dual_key, feature_line])

    header = "frameshift, spliceLength, len_diff, nlenDiff, len_after, coreC1, coreA, coreC2, canCore, disorderC1, " \
             "disorderA, disorderC2, canDisorder, otherDisorder, otherADis, domainC1, domainA, domainC2, canDomain, " \
             "otherDomain, otherADom, enzymatic, elmC1, elmA, elmC2, canElm, mutationC1, mutationA, mutationC2, " \
             "canMutation, ptmC1, ptmA, ptmC2, canPtm, transmemC1, transmemA, transmemC2, canTransmem, degree, " \
             "seqConAve, seqConMin, seqConMax, classifier, spSpec, consAS, gspec, ts"

    print >> write_all, "IE, " + header
    for t in output_positiveI:
        print >> write_all_names, t[0]
        print >> write_all, "1, " + t[1]

    for t in output_positiveE:
        print >> write_all_names, t[0]
        print >> write_all, "0, " + t[1]
