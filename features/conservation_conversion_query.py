# creates query file for converting between hg19 and 18


def create_query_file(as_location_file, conservation_query_output_location):
    as_location = open(as_location_file, 'r')
    write_to = open(conservation_query_output_location, 'w')

    for line in as_location:
        tokens = line.split()
        chromosome = tokens[0]
        start = tokens[1]
        end = tokens[2]
        # asid = tokens[3]
        # type = tokens[4]
        # strand = tokens[5]

        output = "chr" + chromosome + ":" + start + "-" + end
        print >> write_to, output
    write_to.close()
