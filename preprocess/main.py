# Main preprocessing script
import sys
import json
import pickle
from helpers import normalize_unicode_data, load_sequence, load_assembled_transcripts, \
    process_transcripts, fetch_events
from reference_genome import ReferenceGenome

if __name__ == "__main__":

    # Import preprocess settings
    preprocess_settings = json.load(open('preprocess_settings.json', 'r'))

    ##########################################################################

    # # LOAD REFERENCE GENOME
    # ref_genome = ReferenceGenome()
    # REFERENCE_GENOME_DIRECTORY = normalize_unicode_data(
    #     preprocess_settings["REFERENCE_GENOME_DIRECTORY"]
    # )
    #
    # for chromosome in ref_genome.chromosomes:
    #     sequence = load_sequence(REFERENCE_GENOME_DIRECTORY, chromosome)
    #     ref_genome.add_genome(chromosome, sequence)
    #
    # # PICKLE REFERENCE GENOME
    # print "PICKLING REFERENCE GENOME"
    # with open('../input/GRCh37_refGenome_pickle.pkl', 'wb') as output:
    #     pickle.dump(ref_genome, output, pickle.HIGHEST_PROTOCOL)

    # LOAD PICKLED REFERENCE GENOME
    print "LOADING PICKLED REFERENCE GENOME"
    with open('../input/GRCh37_refGenome_pickle.pkl', 'rb') as object_to_load:
        ref_genome = pickle.load(object_to_load)
    print "PICKLED REFERENCE GENOME LOADED"

    ##########################################################################

    # READ TRANSCRIPT TO FIND ALTERNATIVE SPLICING EVENTS

    print "READING TRANSCRIPT FILES..."
    list_of_transcript_files = []
    for transcript_location in preprocess_settings["LIST_OF_TRANSCRIPTS"]:
        list_of_transcript_files.append(normalize_unicode_data(transcript_location))

    dict_of_transcripts = {}

    for transcript_file in list_of_transcript_files:
        dict_of_transcripts[transcript_file] = load_assembled_transcripts(
            transcript_file,
            ref_genome
        )

    ##########################################################################

    # CREATE A DICTIONARY OF ALL UNIQUE TRANSCRIPTS

    FPKM_THRESHOLD = float(normalize_unicode_data(preprocess_settings["FPKM_THRESHOLD"]))

    processed_transcripts = process_transcripts(list_of_transcript_files, dict_of_transcripts,
                                                FPKM_THRESHOLD)
    dictionary_of_unique_transcripts = processed_transcripts["dictionary_of_unique_transcripts"]
    list_transcripts = processed_transcripts["list_transcripts"]

    ##########################################################################

    # group the transcript to find the alternative events

    # I used to group by gene_id, but this become problematic
    # when I take data from different readings (use diff id for
    # the same transcripts when no reference gene was found)
    # then use TSS to group (transcript start splicing)

    dict_group_transcripts = {}

    # for the Novel isoforms, the gene id is not consistent, then I will group
    # isoforms for chromosome and TSS (transcript start site)

    for transcript in dictionary_of_unique_transcripts.itervalues():
        # TSS id
        TSS_id = '-' + transcript.chromosome + '-' + str(transcript.exons[0].start)
        # if not dict_group_transcripts.has_key(transcript.gene_id):
        if TSS_id not in dict_group_transcripts:
            dict_group_transcripts[TSS_id] = [transcript]
        else:
            dict_group_transcripts[TSS_id].append(transcript)

    # extract each transcript of this gene

    aslocation = open('ASlocation.out', 'w')
    complete = open('complete_transcripts.fasta', 'w')

    for TSS_id, list_transcripts in dict_group_transcripts.iteritems():
        if len(list_transcripts) > 1:
            fetch_events(list_transcripts, aslocation, complete)
