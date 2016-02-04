# reads uniprot transmembrane file and transmembrane region features
import os
import sys


#################
#### FUNCTIONS

def score_differences(mapping, uniprot, start, end):
    try:
        count = 0
        if mapping.has_key(uniprot):
            if start <= end:
                for i in range(start, end + 1):
                    if not mapping[uniprot].has_key(i):
                        pass
                    elif mapping[uniprot][i] == "*":
                        count = count + 1
                return (1.0 * count) / (end - start + 1)

    except KeyError, e:
        print 'WARNING!: Features not found '
        print  str(e)
        print 'Please check ID tags or/and generate the missing features'

        return count # count is already referenced before...


#################
#### PATH VARIABLE

try:
    PulsePATH = os.environ['PULSE_PATH']
except KeyError:
    print 'ENVIOROMENT VARIABLE PULSE_PATH not set up'
    print 'TRYNG WITH WORKING DIR'
    PulsePATH = os.getcwd()

#################
##### INPUT FILES
try:
    # Load Transmembrane info
    uniprot_TM_indices_DB = open(PulsePATH + '/param_files/F_uniprot_transmem_indices.inp', 'r')

    # Read the file with the mapping between a splicing and its anchor
    map_file = sys.argv[1]
    readFrom1 = open(map_file, 'r')

except IndexError:
    print 'Map file must be provided'
    sys.exit()


###############
#### OUTPUT FILES

writeTo = open('transmem_read.out', 'w')

###############
#### Init VARs

uniprot_to_index_to_whatever = {}


###############
#### MAIN


for line in uniprot_TM_indices_DB:
    tokens = line.split("\t")

    try:
        uniprot = tokens[0].strip()
        start = int(tokens[2].strip())
        end = int(tokens[3].strip())

        for i in range(start, end + 1):

            if uniprot_to_index_to_whatever.has_key(uniprot):
                uniprot_to_index_to_whatever[uniprot][i] = "*"
            else:
                uniprot_to_index_to_whatever[uniprot] = {i: "*"}

    except ValueError:
        print "Cannot parse: " + line[0:len(line) - 1]

for line in readFrom1:
    tokens = line.split()
    # PARSING ID
    # Customize the line below if is need it
    # it depends of how it is the label of th events
    # This is ready for the default format label splicing plus lenght of the exons -> X.LABEL-###_33=33=33
    asid = tokens[0].split("_")[0]
    prot = tokens[1]
    # PARSING ID
    sstart = int(tokens[2])
    start = int(tokens[3])
    end = int(tokens[4])
    eend = int(tokens[5])

    c1Count = 0
    aCount = 0
    c2Count = 0
    canonicalAbsolute = 0

    if not uniprot_to_index_to_whatever.has_key(prot):
        c1Count = 0
        aCount = 0
        c2Count = 0
        canonicalAbsolute = 0
        otherAbsolute = 0
    else:
        c1Count = score_differences(uniprot_to_index_to_whatever, prot, sstart, start)
        aCount = score_differences(uniprot_to_index_to_whatever, prot, start, end)
        c2Count = score_differences(uniprot_to_index_to_whatever, prot, end, eend)

        # use teh new prot length provided
        protLen = int(line.split("\t")[7].strip())

        canonicalAbsolute = score_differences(uniprot_to_index_to_whatever, prot, 1, protLen)

    print >> writeTo, tokens[0] + "\t" + prot + "\t" + repr(c1Count) + "\t" + repr(aCount) + "\t" + repr(
        c2Count) + "\t" + repr(canonicalAbsolute)  # + "\t" + repr(otherAbsolute)
