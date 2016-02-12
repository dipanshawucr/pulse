# generates sequence conservation features
import sys
import os

#################
# PATH VARIABLE

try:
    PulsePATH = os.environ['PULSE_PATH']
except KeyError:
    print 'ENVIOROMENT VARIABLE PULSE_PATH not set up'
    print 'TRYNG WITH WORKING DIR'
    PulsePATH = os.getcwd()

#################
# INPUT FILES

try:
    readFrom2 = open(PulsePATH + '/src_features/param_files/F_phastCons.hg18.bed.inp', 'r')
    as_location_file = sys.argv[1]
    remapped_coordinates = sys.argv[2]
    as_location_file_object = open(as_location_file, 'r')
    remapped_coordinates_object = open(remapped_coordinates, 'r')

except IndexError:
    print 'Alternative Splicing Location and the Remap hg18 file must be provided'
    sys.exit()

#################
# OUPUT FILES

writeTo = open('sequenceCon_read.out', 'w')

#################
# Init VARs

old_to_new = {}
new_to_score = {}

#################
# MAIN

for line in remapped_coordinates_object:
    print line
    if line[0] != '#':
        tokens = line.split(',')
        old_key = tokens[3] + "//" + tokens[7] + "//" + tokens[8]
        new_key = tokens[4] + "//" + repr(int(tokens[12]) - 1) + "//" + tokens[13]
        old_to_new[old_key] = new_key

for line in readFrom2:
    if line[0] != '#':
        tokens = line.split()
        new_key = tokens[0] + "//" + tokens[1] + "//" + tokens[2]
        new_to_score[new_key] = tokens[6] + "\t" + tokens[7] + "\t" + tokens[8]

for line in as_location_file_object:
    tokens = line.split()

    chromosome = tokens[0]
    start = tokens[1]
    end = tokens[2]
    asid = tokens[3]
    type = tokens[4]
    strand = tokens[5]

    try:
        key = "chr" + chromosome + "//" + start + "//" + end
        new_key = old_to_new[key]
        score = new_to_score[new_key]

        if type == "A":
            print >> writeTo, asid + "\t" + score

    except KeyError:
        print key + "\t" + new_key
