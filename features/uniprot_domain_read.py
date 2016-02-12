# reads uniprot domain file and generates domain features
from features_helpers import link_db, score_differences_pfam
import time


def fetch_feature(feature, anchor, path_ddbb):
    (c, db) = link_db(path_ddbb)
    pfam_hits = []
    query = c.execute(""" SELECT * FROM %s WHERE id='%s' """ % (feature, anchor))
    # query = c.execute("""SELECT * FROM uniprot_ptms WHERE uniprot_id = "Q8IZP0" """)
    # Row_id Uniprot_id postion Value ['*'] in pfam ['*',enzymatic bool]
    for row in query:
        uniprot_id = row[0]
        start = int(row[1])
        end = int(row[2])
        pfam_id = row[3]
        pfam_hits.append([uniprot_id, start, end, pfam_id])
    c.close()
    return pfam_hits


def extract_from_pfam_special_db(f_pfam_special_db):
    pfam_special = {}
    for line in f_pfam_special_db:
        tokens = line.split()
        pfam = tokens[1]
        pfam_special[pfam] = True
    f_pfam_special_db.close()
    return pfam_special


def build_uniprot_to_index_domain(pfamscan_output_isoforms, pfam_special):
    uniprot_to_index_to_domain = {}
    for line in pfamscan_output_isoforms:
        if line[0] != "#" and line.strip() != "":
            tokens = line.split()
            try:
                asid = tokens[0]
                start = int(tokens[1])
                end = int(tokens[2])
                pfam = tokens[5].split(".")[0]
                enzymatic = pfam in pfam_special
                for i in range(start, end + 1):
                    if asid in uniprot_to_index_to_domain:
                        uniprot_to_index_to_domain[asid][i] = ["*", enzymatic]
                    else:
                        uniprot_to_index_to_domain[asid] = {i: ["*", enzymatic]}
            except ValueError:
                print "Cannot parse: " + line[0:len(line) - 1]
            except IndexError:
                print "Index Error: " + line.strip()
    return uniprot_to_index_to_domain


def get_uniprot_domain_read(f_pfam_special_db_location, canonical_db_location, map_file, pfamscan_output,
                            domain_read_output_location):
    f_pfam_special_db = open(f_pfam_special_db_location, 'r')
    pfamscan_output_isoforms = open(pfamscan_output, 'r')
    write_to = open(domain_read_output_location, 'w')

    pfam_special = extract_from_pfam_special_db(f_pfam_special_db)
    uniprot_to_index_to_domain = build_uniprot_to_index_domain(pfamscan_output_isoforms, pfam_special)

    index_file = open(map_file, 'r')
    for line in index_file:
        tokens = line.split()
        asid = tokens[0]
        prot = tokens[1]
        pfams = fetch_feature('uniprot_pfams', prot, canonical_db_location)
        for hit in pfams:
            start = hit[1]
            end = hit[2]
            pfam = hit[3]
            enzymatic = pfam in pfam_special
            for i in range(start, end + 1):
                if prot in uniprot_to_index_to_domain:
                    uniprot_to_index_to_domain[prot][i] = ["*", enzymatic]
                else:
                    uniprot_to_index_to_domain[prot] = {i: ["*", enzymatic]}
    index_file.close()

    # Generate Features
    index_file = open(map_file, 'r')
    for line in index_file:
        tokens = line.split()
        asid = tokens[0]
        prot = tokens[1]
        sstart = int(tokens[2])
        start = int(tokens[3])
        end = int(tokens[4])
        eend = int(tokens[5])

        rough_a_length = int(int(tokens[0].split("_")[-1].split("=")[1]) / 3)

        if asid[0] == "I":
            rough_a_length = 0

        if prot not in uniprot_to_index_to_domain:
            c1_count = 0
            a_count = 0
            c2_count = 0
            canonical_absolute = 0
            other_absolute = 0
            other_a = 0
            enzymatic = 0
        else:
            c1_count = score_differences_pfam(uniprot_to_index_to_domain, prot, sstart, start)[0]
            a_count = score_differences_pfam(uniprot_to_index_to_domain, prot, start, end)[0]
            c2_count = score_differences_pfam(uniprot_to_index_to_domain, prot, end, eend)[0]
            enzymatic = score_differences_pfam(uniprot_to_index_to_domain, prot, start, end)[1]

            prot_len = int(line.split("\t")[7].strip())
            other_prot = asid
            other_prot_len = int(line.split("\t")[8].strip())
            canonical_absolute = score_differences_pfam(uniprot_to_index_to_domain, prot, 1, prot_len)[0]
            other_absolute = score_differences_pfam(uniprot_to_index_to_domain, other_prot, 1, other_prot_len)[0]

            other_a_end = start + rough_a_length
            if other_a_end > other_prot_len:
                other_a_end = other_prot_len

            other_a = score_differences_pfam(uniprot_to_index_to_domain, other_prot, start, other_a_end)[0]

        print >> write_to, tokens[0] + "\t" + prot + "\t" + repr(c1_count) + "\t" + repr(a_count) + "\t" + repr(
            c2_count) + "\t" + repr(canonical_absolute) + "\t" + repr(other_absolute) + "\t" + repr(other_a) + "\t" + \
                           repr(enzymatic)

    write_to.close()
