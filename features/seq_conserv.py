# generates sequence conservation features


def generate_old_to_new(remapped_coordinates_object):
    old_to_new = {}
    for line in remapped_coordinates_object:
        if line[0] != '#':
            tokens = line.split('\t')
            old_key = tokens[3] + "//" + tokens[7] + "//" + tokens[8]
            new_key = tokens[4] + "//" + repr(int(tokens[12]) - 1) + "//" + tokens[13]
            old_to_new[old_key] = new_key
    return old_to_new


def generate_new_to_score(f_phastcons_db):
    new_to_score = {}
    for line in f_phastcons_db:
        if line[0] != '#':
            tokens = line.split()
            new_key = tokens[0] + "//" + tokens[1] + "//" + tokens[2]
            new_to_score[new_key] = tokens[6] + "\t" + tokens[7] + "\t" + tokens[8]
    return new_to_score


def generate_sequence_conservation_features(f_phastcons_db_location, as_location_file,
                                            remapped_coordinates, output_location):
    f_phastcons_db = open(f_phastcons_db_location, 'r')
    as_location_file_object = open(as_location_file, 'r')
    remapped_coordinates_object = open(remapped_coordinates, 'r')

    write_to = open(output_location, 'w')

    old_to_new = generate_old_to_new(remapped_coordinates_object)
    new_to_score = generate_new_to_score(f_phastcons_db)

    for line in as_location_file_object:
        tokens = line.split()
        chromosome = tokens[0]
        start = tokens[1]
        end = tokens[2]
        asid = tokens[3]
        type = tokens[4]
        # strand = tokens[5]
        try:
            key = "chr" + chromosome + "//" + start + "//" + end
            new_key = old_to_new[key]
            score = new_to_score[new_key]
            if type == "A":
                print >> write_to, asid + "\t" + score
        except KeyError:
            print key + "\t" + new_key
