# generates event conservation feature table
from seq_conserv import generate_old_to_new


def conv(value):
    if value:
        return "1"
    else:
        return "0"


def generate_new_to_nm(f_phastcons_db):
    new_to_nm = {}
    for line in f_phastcons_db:
        tokens = line.split()
        new_key = tokens[0] + "//" + tokens[1] + "//" + tokens[2]
        new_to_nm[new_key] = tokens[3]
    return new_to_nm


def generate_nm_to_quad(f_event_conservation_db):
    nm_to_quad = {}
    for line in f_event_conservation_db:
        tokens = line.split()
        key = tokens[0]
        quad = map(conv, ["Classifier" in tokens[1], "Sp-spec" in tokens[1], "Cons-AS" in tokens[1], "G-spec" in tokens[1]])
        nm_to_quad[key] = quad
    return nm_to_quad


def generate_event_conservation_feature_table(f_phastcons_db_location, f_event_conservation_db_location,
                                              as_location_file, remapped_coordinates_location,
                                              event_conservation_output):
    f_phastcons_db = open(f_phastcons_db_location, 'r')
    f_event_conservation_db = open(f_event_conservation_db_location, 'r')
    as_location_obj = open(as_location_file, 'r')
    remapped_coordinates_obj = open(remapped_coordinates_location, 'r')

    write_to = open(event_conservation_output, 'w')

    old_to_new = generate_old_to_new(remapped_coordinates_obj)
    new_to_nm = generate_new_to_nm(f_phastcons_db)
    nm_to_quad = generate_nm_to_quad(f_event_conservation_db)

    for line in as_location_obj:
        tokens = line.split()
        chromosome = tokens[0]
        start = tokens[1]
        end = tokens[2]
        asid = tokens[3]
        type = tokens[4]
        # strand = tokens[5]
        if type == "A":
            try:
                key = "chr" + chromosome + "//" + start + "//" + end
                new_key = old_to_new[key]
                nm = new_to_nm[new_key]
                quad = ["0", "0", "0", "0"]
                if nm in nm_to_quad:
                    quad = nm_to_quad[nm]
                print >> write_to, asid + "\t" + quad[0] + "\t" + quad[1] + "\t" + quad[2] + "\t" + quad[3]
            except KeyError:
                print key + "\t" + new_key
