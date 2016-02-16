# reads isoforms index, return a file with all that need salble calcs


try:
    map_file = sys.argv[1]
    readFrom1 = open(map_file, 'r')
    insable_file = PulsePATH+'/src_features/param_files/insable.list'
    pfam_file = PulsePATH+'/src_features/param_files/pfam_scan_uniprot.inp'
except IndexError:
    print 'Map and iupred output files must be provided'
    sys.exit()

dict_canonical = {}
filename = '/home/ccorbi/Work/ddbb_fasta/Uniprot/uniprot_sprot.fasta'
with open(filename, 'r') as input_file:
    for line in input_file:

        if line[0] == '>':
            protein_id = line.split('|')[1]

        else:
            if protein_id in dict_canonical:
                dict_canonical[protein_id] += line.strip()
            else:
                dict_canonical[protein_id] = line.strip()


#################
# OUTPUT FILES

writesable = open('need_sable.out', 'w')
writepfam = open('need_pfam.out', 'w')

# Init VARs

uniprot_to_index_to_disorder = {}
insable = {}
inpfam = {}

needs = {}
needp = {}

###############
# MAIN

with open(insable_file, 'r') as input_file:
    for line in input_file:
        data = line.strip()
        insable[data] = 1

with open(pfam_file, 'r') as input_file:
    for line in input_file:
        data = line.strip('\t')
        inpfam[data[0]] = 1


print len(dict_canonical)
for line in readFrom1:
    tokens = line.split('\t')
    asid = tokens[0].split("_")[0]
    prot = tokens[1]
    # print prot
    if prot in insable:
        pass
    else:
        if prot in dict_canonical:
            if prot not in needs:
                print >> writesable, '>' + prot
                print >> writesable, dict_canonical[prot]
                needs[prot] = 1
        else:
            print '>'+prot

    if prot in inpfam:
        pass
    else:
        if prot in dict_canonical:
            if prot not in needp:
                print >> writepfam, '>' + prot
                print >> writepfam, dict_canonical[prot]
                needp[prot] = 1
        else:
            print '>' + prot
