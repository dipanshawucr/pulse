# applies length filter on the already gene filtered query results
import sys
import os


def pbegend(pstart, pend, plength):
    phend = (pend + threshold / 3 > plength)
    phbeg = (pstart - threshold / 3 < 1)

    return [phbeg, phend]


def mbegend(nstart, nend, nlength, backwards):
    if backwards:
        temp = nstart
        nstart = nend
        nend = temp
    mhend = (nend + threshold > nlength)
    mhbeg = (nstart - threshold < 1)

    return [mhbeg, mhend]


#################
##### INPUT FILES

try:
    blast_file = sys.argv[1]
    readFrom = open(blast_file, 'r')
    uniprot_fasta = sys.argv[2]
    isoform_fasta = sys.argv[3]
except:
    print 'Blast output file, uniprot in fasta and isoform seq in fasta must be provided'
    sys.ext()

try:
    os.makedirs("temp")
except OSError, e:
    if e.errno != 17:
        raise
    print e
    pass

writeTo = open('./temp/rnaseq_huniprot_corrected_len.txt', 'w')

threshold = 9  # nucleotides
pidentThreshold = 95.0
ignoreExclusion = False
dealWithNucleotide = False
lineNumber = 0

for line in readFrom:

    (queryId, subjectId, percIdentity, alnLength, mismatchCount, gapOpenCount, queryStart, queryEnd, subjectStart,
     subjectEnd, eVal, bitScore) = line.split('\t')
    try:
        # PARSING ID
        # qId = queryId.split("-")[-1]
        # print str(queryId.split("-")[-1:]),type(queryId.split("-")[1])
        # PARSING ID
        snucleotideLengths = queryId.split("-")[-1].split("=")
        nucleotideLengths = [int(snucleotideLengths[0]), int(snucleotideLengths[1]), int(snucleotideLengths[2])]
    # print nucleotideLengths
    # print snucleotideLengths
    except:
        print 'WARNING WITH' + queryId

    pstart = int(subjectStart)  # /3.0
    pend = int(subjectEnd)  # /3.0
    plength = int(alnLength)  # /3.0

    if dealWithNucleotide:
        pstart /= 3.0
        pend /= 3.0
        plength /= 3.0

    nstart = int(queryStart)
    nend = int(queryEnd)
    nlength = sum(nucleotideLengths)

    pident = float(percIdentity)

    if queryId[0] == "E":
        if ignoreExclusion:
            continue
        nlength -= nucleotideLengths[1]

    toPrint = True

    backwards = (nstart > nend)

    [hitPBeg, hitPEnd] = pbegend(pstart, pend, plength)
    [hitNBeg, hitNEnd] = mbegend(nstart, nend, nlength, backwards)

    existBeg = True
    existEnd = True

    if not hitNBeg and not hitPBeg:
        toPrint = False
        existBeg = False

    if not hitNEnd and not hitPEnd:
        toPrint = False
        existEnd = False

    if not (existBeg or existEnd):
        print line.strip() + " no_both"
    elif not existBeg:
        if backwards:
            print line.strip() + " no_c2"
        else:
            print line.strip() + " no_c1"
    elif not existEnd:
        if backwards:
            print line.strip() + " no_c1"
        else:
            print line.strip() + " no_c2"
    else:
        hitLength = nend - nstart + 1
        if hitLength < 0:
            hitLength *= -1
            hitLength += 2

        if ((hitLength + threshold) / 3 < pend - pstart + 1) or ((hitLength - threshold) / 3 > pend - pstart + 1):
            toPrint = False
            print line.strip() + " missing_chunks"

    if pident < pidentThreshold:
        toPrint = False
        print line.strip() + " pident_failure"

    if queryId[0] == "E":
        # make sure A site is within
        if not backwards:
            if not (nucleotideLengths[0] <= nend and nucleotideLengths[0] >= nstart):
                toPrint = False
            # print "SHANKED: " + line.strip()
        else:
            if not (nucleotideLengths[0] <= nstart and nucleotideLengths[0] >= nend):
                # print "SHANKED: " + line.strip()
                toPrint = False

    if toPrint:
        print >> writeTo, line.strip()[0:len(line)]

# filters multiple uniprot hits by taking the best one in terms of evalue

readFrom = open('./temp/rnaseq_huniprot_corrected_len.txt', 'r')
writeTo = open('./temp/rnaseq_huniprot_corrected_len_collapsed.txt', 'w')

oldID = ""

for line in readFrom:
    tokens = line.split('\t')
    # PARSING ID
    asid = tokens[0]
    prot = tokens[1]
    # PARSING ID
    # id = asid + "//" + prot
    id = asid

    if oldID != id:
        print >> writeTo, line.strip()

    oldID = id

###
# LOAD UNIPROT
uniprot_ddbb = {}
seq = ''
uid = '|'
SEP = '|'
# print uniprot_fasta
with open(uniprot_fasta, 'r') as input_file:
    for line in input_file:
        if line[0] == '>':
            uniprot_ddbb[uid] = len(seq)
            # PARSING ID
            uid = line.split(SEP)[1]
            seq = ''
        # PARSING ID
        else:
            seq += line.strip()

uniprot_ddbb[uid] = len(seq)
###
# LOAD ISOFORM IN FASTA
isoforms_ddbb = {}
seq = ''
uid = ''

SEP = '>'
with open(isoform_fasta, 'r') as input_file:
    for line in input_file:
        if line[0] == '>':
            isoforms_ddbb[uid] = len(seq)
            # CHANGES THE LINE BELOW if is need it
            # PARSING ID
            uid = line.split(SEP)[1].strip()
            seq = ''
        # PARSING ID
        else:
            seq += line.strip()

isoforms_ddbb[uid] = len(seq)

# tags on spliced exon indices on query hits using pstart pend nstart and nend along with A exon positions
# this information is obtained from the splice sequence to uniprot mapping file


readFrom = open('./temp/rnaseq_huniprot_corrected_len_collapsed.txt', 'r')
writeTo = open('uniprot_exon_indices_map.out', 'w')

hasMapping = {}

for line in readFrom:
    tokens = line.split('\t')
    # PARSING ID
    id = tokens[0]  # asid
    uniprot = tokens[1].split("|")[1]
    # PARSING ID
    if hasMapping.has_key(id):
        hasMapping[id].append(uniprot)
    else:
        hasMapping[id] = [uniprot]

readFrom.close()

readFrom = open('./temp/rnaseq_huniprot_corrected_len_collapsed.txt', 'r')

for line in readFrom:
    tokens = line.split('\t')

    # PARSING ID
    id = tokens[0]  # asid
    uniprot = tokens[1].split("|")[1]
    # PARSING ID

    snucleotideLengths = tokens[0].split("-")[-1].split("=")
    nucleotideLengths = [int(snucleotideLengths[0]), int(snucleotideLengths[1]), int(snucleotideLengths[2])]

    pstart = int(tokens[8])
    pend = int(tokens[9])
    nstart = int(tokens[6])
    nend = int(tokens[7])

    if nstart < nend:  # these could be negative
        c1HitLength = nucleotideLengths[0] - nstart + 1
        c2HitLength = nucleotideLengths[2] - (sum(nucleotideLengths) - nend + 1)
    else:
        c1HitLength = nucleotideLengths[0] - nend + 1
        c2HitLength = nucleotideLengths[2] - (sum(nucleotideLengths) - nstart + 1)

    if c1HitLength < 0:
        c1HitLength = 0;
    if c2HitLength < 0:
        c2HitLength = 0;

    aStart = pstart + c1HitLength / 3
    aEnd = pend - c2HitLength / 3
    # PARSING ID
    # ID HANDELING change it if you need it
    transcript_isoform = id[2:].split('-')[0]
    # PARSING ID

    try:

        if aStart > 0 and aEnd > 0 and aEnd >= aStart and "-" not in uniprot:
            if id[0] == "I":
                listProt = []
                if hasMapping.has_key("E" + id[1:len(id)]):
                    listProt = hasMapping["E" + id[1:len(id)]]
                print >> writeTo, id + "\t" + uniprot + "\t" + repr(pstart) + "\t" + repr(aStart) + "\t" + repr(
                    aEnd) + "\t" + repr(pend) + "\t" + repr(listProt) + "\t" + str(uniprot_ddbb[uniprot]) + "\t" + str(
                    isoforms_ddbb[transcript_isoform])
            else:
                listProt = []
                if hasMapping.has_key("I" + id[1:len(id)]):
                    listProt = hasMapping["I" + id[1:len(id)]]
                print >> writeTo, id + "\t" + uniprot + "\t" + repr(pstart) + "\t" + repr(aStart) + "\t" + repr(
                    aStart) + "\t" + repr(pend) + "\t" + repr(listProt) + "\t" + str(
                    uniprot_ddbb[uniprot]) + "\t" + str(isoforms_ddbb[transcript_isoform])
        else:
            pass

    except KeyError, e:
        print '%s NOT FOUND ' % str(e)
        print 'Please check your Uniprot file or/and your poll of isoforms'
    # print id
