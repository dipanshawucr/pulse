# generates network features using gene names


def generate_dict_uniprot(file_gene_id_index):
    # LOAD TABLE Gene Name - Uniprot_ID
    # Return a Dict key:Uniprot Value: GeneWiki Name
    dict_gene_wiki = {}
    for line in file_gene_id_index:
        tab_data = line.split('\t')
        dict_gene_wiki[tab_data[1].strip()] = tab_data[0]
    return dict_gene_wiki


def generate_network_features(f_uniprot_genewiki_location, f_degree_location, map_file, output_location):
    file_gene_id_index = open(f_uniprot_genewiki_location, 'r')
    f_degree_db = open(f_degree_location, 'r')
    file_anchor_map = open(map_file, 'r')
    write_to = open(output_location, 'w')

    dict_uniprot = {}
    name_to_unigene = {}
    dict_gene_wiki = generate_dict_uniprot(file_gene_id_index)

    # LOAD TABLE ISOFORM-ANCHOR_Uniprot
    # Return a Dict key:Uniprot Value:Alt Splt ID
    for line in file_anchor_map:
        tab_data = line.split('\t')
        # Key:Uniprot-Anchor Value:Alt Splic Event ID
        # PARSING ID
        as_id = tab_data[0]
        dict_uniprot[tab_data[1].strip()] = as_id
        # PARSING ID
        # look if there is a uniprot->gene
        for gene_id, uniprot in dict_gene_wiki.iteritems():
            if tab_data[1] == uniprot:
                name_to_unigene[gene_id] = as_id
                break

    # LOAD degree File select matches and print score
    for line in f_degree_db:
        print line
        tokens = line.split(",")

        name = tokens[1].strip()
        score = tokens[2].strip()
        print name, score
        if name in name_to_unigene:
            print "NAME: ", name
            as_id = name_to_unigene[name]
            print >> write_to, as_id + "\t" + score
