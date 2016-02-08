# reads uniprot disorder file and generates disorder features
# do absolute disorder as well
import sys
import os
import sqlite3


def score_differences(mapping, uniprot, start, end):
    try:
        count = 0
        if start <= end:
            for i in range(start, end + 1):
                if mapping[uniprot.strip()][i] == "*":
                    count = count + 1
            count = (1.0 * count) / (end - start + 1)

        return count

    except KeyError, e:
        print 'WARNING!: Features not found '
        print  str(e), uniprot
        print 'Please check ID tags or/and generate the missing features'


def link_db(db_path):
    """
    Initializes database connection

    :param db_path:
    :return:
    """

    try:
        db = sqlite3.connect(db_path)
    except sqlite3.Error, errmsg:
        print 'DB not available ' + str(errmsg)
        sys.exit()

    db_cursor = db.cursor()
    return db_cursor, db


def fetch_feature(feature, anchors_map, path_ddbb):
    # connect ddbb
    # global DDBB
    # print path_ddbb
    (c, db) = link_db(path_ddbb)

    # Init Var
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

if __name__ == "__main__":

    # PATH VARIABLE
    try:
        PulsePATH = os.environ['PULSE_PATH']
    except KeyError:
        print 'ENVIOROMENT VARIABLE PULSE_PATH not set up'
        print 'TRYNG WITH WORKING DIR'
        PulsePATH = os.getcwd()

    # INPUT FILES
    try:
        map_file = sys.argv[1]
        iupred_output = sys.argv[2]
        # readFrom1 = open(map_file, 'r')
        assert open(map_file, 'r')
        assert open(iupred_output, 'r')
        canonical_ddbb = PulsePATH + '/src_features/database/info_canonical_v3.ddbb'
        print "Directory of database being used: " + canonical_ddbb

    except IndexError:
        print 'Map and iupred output files must be provided'
        sys.exit()

    # OUTPUT FILES
    writeTo = open('./output/disorder_read.out', 'w')

    # Init VARs
    uniprot_to_index_to_disorder = {}
    anchors_map = []

    # Load IUPred info from anchor in ddbb
    # Load anchor ids
    readFrom1 = open(map_file, 'r')
    for line in readFrom1:
        tokens = line.split('\t')
        prot = tokens[1]
        anchors_map.append(prot)
    readFrom1.close()

    uniprot_to_index_to_disorder = fetch_feature('uniprot_iupred', anchors_map, canonical_ddbb)

    readFrom = open(iupred_output, 'r')
    for line in readFrom:
        print line
        tokens = line.split('\t')

        if len(tokens) > 1:
            try:
                # PARSING ID
                prot = tokens[0]
                # Customize here if you want to mod the label
                # if "TCONS" in prot:
                # 	prot = prot.split('_')[0]
                # PARSING ID
                index = int(tokens[1])
                disordered = tokens[2].strip()

                if uniprot_to_index_to_disorder.has_key(prot):
                    uniprot_to_index_to_disorder[prot][index] = disordered
                else:
                    uniprot_to_index_to_disorder[prot] = {index: disordered}
            except ValueError:
                print "Cannot parse: " + line[0:len(line) - 1]

    readFrom.close()

    # print uniprot_to_index_to_disorder
    # sys.exit()
    readFrom1 = open(map_file, 'r')
    for line in readFrom1:
        tokens = line.split('\t')
        # print "inside"
        # print >>writeTo,'works'
        # PARSING ID
        # Customize the line below if is need it
        # it depends of how it is the label of th events
        # This is ready for the default format label splicing plus lenght of the exons -> X.LABEL-###_33=33=33
        # PARSING ID
        asid = tokens[0]  # .split("_")[0]
        prot = tokens[1]
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

        if not uniprot_to_index_to_disorder.has_key(prot):
            print prot
            c1Count = 0
            aCount = 0
            c2Count = 0
            canonicalAbsolute = 0
            otherAbsolute = 0
            otherA = 0
        else:
            c1Count = score_differences(uniprot_to_index_to_disorder, prot, sstart, start)
            aCount = score_differences(uniprot_to_index_to_disorder, prot, start, end)
            c2Count = score_differences(uniprot_to_index_to_disorder, prot, end, eend)

            # use teh new prot length provided
            protLen = int(line.split("\t")[7].strip())
            otherProt = asid  # or token[0] Customize this line if is need it for test [asid+'//'+prot+'-ST']
            otherProtLen = int(line.split("\t")[8].strip())

            canonicalAbsolute = score_differences(uniprot_to_index_to_disorder, prot, 1, protLen)
            otherAbsolute = score_differences(uniprot_to_index_to_disorder, otherProt, 1, otherProtLen)

            otherAEnd = start + roughALength
            if otherAEnd > otherProtLen:
                otherAEnd = otherProtLen
            # print line.strip()

            otherA = score_differences(uniprot_to_index_to_disorder, otherProt, start, otherAEnd)

        print >> writeTo, tokens[0] + "\t" + prot + "\t" + repr(c1Count) + "\t" + repr(aCount) + "\t" + repr(
            c2Count) + "\t" + repr(canonicalAbsolute) + "\t" + repr(otherAbsolute) + "\t" + repr(otherA)
