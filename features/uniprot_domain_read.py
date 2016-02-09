# reads uniprot domain file and generates domain features
import sys
from features_helpers import link_db, score_differences_pfam


def fetch_feature(feature, anchor, path_ddbb):
    # connect ddbb
    # global DDBB
    # print path_ddbb
    (c, db) = link_db(path_ddbb)

    # Init Var
    pfam_hits = []
    # for anchor in anchors_map:
    # print anchor
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


def get_uniprot_domain_read(f_pfam_special_db_location, canonical_db_location):

    try:
        f_pfam_special_db = open(f_pfam_special_db_location, 'r')
        # readFrom3 = open(PulsePATH+'/param_files/pfam_scan_uniprot.inp', 'r')
        map_file = sys.argv[1]
        pfamscan_output = sys.argv[2]

        assert open(map_file, 'r')
        pfamscan_output_isoforms = open(pfamscan_output, 'r')
        canonical_ddbb = PulsePATH + '/src_features/database/info_canonical_v3.ddbb'
    # canonical_ddbb = PulsePATH+'/database/pulse_canonical_v2.ddbb'

    except IndexError:
        print 'Map and pfamscan output files must be provided'
        sys.exit()

    # OUTPUT FILES
    writeTo = open('domain_read.out', 'w')

    # Init VARs
    uniprot_to_index_to_domain = {}
    pfam_special = {}

    # MAIN

    for line in f_pfam_special_db:
        tokens = line.split()
        pfam = tokens[1]
        pfam_special[pfam] = True
    f_pfam_special_db.close()

    for line in pfamscan_output_isoforms:
        tokens = line.split()
        # print tokens
        try:
            asid = tokens[0]

            # if "|" in uniprot:
            # 	uniprot = uniprot.split("|")[1]

            start = int(tokens[1])
            end = int(tokens[2])

            pfam = tokens[5].split(".")[0]
            enzymatic = pfam_special.has_key(pfam)

            for i in range(start, end + 1):
                if uniprot_to_index_to_domain.has_key(asid):
                    uniprot_to_index_to_domain[asid][i] = ["*", enzymatic]
                else:
                    uniprot_to_index_to_domain[asid] = {i: ["*", enzymatic]}

        except ValueError:
            print "Cannot parse: " + line[0:len(line) - 1]
        except IndexError:
            print "Index Error: " + line.strip()

    print len(uniprot_to_index_to_domain)
    # print 'normal'
    # load pfam info anchors
    index_file = open(map_file, 'r')
    for line in index_file:
        tokens = line.split()
        # PARSING ID
        # Customize the line below if is need it
        # it depends of how it is the label of th events
        # This is ready for the default format label splicing plus lenght of the exons -> X.LABEL-###_33=33=33

        asid = tokens[0]  # .split("_")[0]
        prot = tokens[1]
        pfams = fetch_feature('uniprot_pfams', prot, canonical_ddbb)
        for hit in pfams:
            start = hit[1]
            end = hit[2]
            pfam = hit[3]
            enzymatic = pfam_special.has_key(pfam)

            for i in range(start, end + 1):
                if uniprot_to_index_to_domain.has_key(prot):
                    uniprot_to_index_to_domain[prot][i] = ["*", enzymatic]
                else:
                    uniprot_to_index_to_domain[prot] = {i: ["*", enzymatic]}

    index_file.close()
    print len(uniprot_to_index_to_domain)
    # Generate Features
    index_file = open(map_file, 'r')
    for line in index_file:
        tokens = line.split()
        # PARSING ID
        # Customize the line below if is need it
        # it depends of how it is the label of th events
        # This is ready for the default format label splicing plus lenght of the exons -> X.LABEL-###_33=33=33

        asid = tokens[0]  # .split("_")[0]
        prot = tokens[1]
        # PARSING ID
        sstart = int(tokens[2])
        start = int(tokens[3])
        end = int(tokens[4])
        eend = int(tokens[5])

        roughALength = int(int(tokens[0].split("_")[-1].split("=")[1]) / 3)

        if asid[0] == "I":
            roughALength = 0

        c1Count = 0
        aCount = 0
        c2Count = 0
        canonicalAbsolute = 0
        otherAbsolute = 0
        otherA = 0

        if not uniprot_to_index_to_domain.has_key(prot):
            c1Count = 0
            aCount = 0
            c2Count = 0
            canonicalAbsolute = 0
            otherAbsolute = 0
            otherA = 0
            enzymatic = 0
        else:
            c1Count = score_differences_pfam(uniprot_to_index_to_domain, prot, sstart, start)[0]
            aCount = score_differences_pfam(uniprot_to_index_to_domain, prot, start, end)[0]
            c2Count = score_differences_pfam(uniprot_to_index_to_domain, prot, end, eend)[0]

            enzymatic = score_differences_pfam(uniprot_to_index_to_domain, prot, start, end)[1]

            # use teh new prot length provided
            protLen = int(line.split("\t")[7].strip())
            otherProt = asid  # or token[0] Customize this line if is need it for test [asid+'//'+prot+'-ST']
            otherProtLen = int(line.split("\t")[8].strip())
            canonicalAbsolute = score_differences_pfam(uniprot_to_index_to_domain, prot, 1, protLen)[0]
            otherAbsolute = score_differences_pfam(uniprot_to_index_to_domain, otherProt, 1, otherProtLen)[0]

            otherAEnd = start + roughALength
            if otherAEnd > otherProtLen:
                otherAEnd = otherProtLen

            otherA = score_differences_pfam(uniprot_to_index_to_domain, otherProt, start, otherAEnd)[0]

        print >> writeTo, tokens[0] + "\t" + prot + "\t" + repr(c1Count) + "\t" + repr(aCount) + "\t" + repr(
            c2Count) + "\t" + repr(canonicalAbsolute) + "\t" + repr(otherAbsolute) + "\t" + repr(otherA) + "\t" + repr(
            enzymatic)
