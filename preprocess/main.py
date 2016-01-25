# Main preprocessing script
import sys
import json
import pickle
from helpers import normalize_unicode_data, load_sequence, load_assembled_transcripts
from reference_genome import ReferenceGenome

if __name__ == "__main__":

    # Import preprocess settings
    preprocess_settings = json.load(open('preprocess_settings.json', 'r'))

    ##########################################################################

    # # LOAD REFERENCE GENOME
    # ref_genome = ReferenceGenome()
    # reference_genome_directory = normalize_unicode_data(
    #     preprocess_settings["REFERENCE_GENOME_DIRECTORY"]
    # )
    #
    # for chromosome in ref_genome.chromosomes:
    #     sequence = load_sequence(reference_genome_directory, chromosome)
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

    print dict_of_transcripts

    ##########################################################################

    dict_uniq_transcrip = {}
    FPKM_THRESHOLD = 0.0

    # for transcript_file in list_transcript_files:
    #     list_transcripts = dict_transcripts[transcript_file]
    #     for transcript in list_transcripts:
    #         exon_ids = ''
    #         for exon in transcript.exons:
    #             exon_ids += str(exon.start) + '-' + str(exon.end) + '.'
    #
    #         transcript_unique_id = transcript.chromosome + '-' + exon_ids
    #         # transcript id is unique for this datasets NO! this is not!! gene_id and transcription_id, not usefull
    #         # transcript_unique_id = transcript.id
    #
    #         if not dict_uniq_transcrip.has_key(transcript_unique_id) and transcript.fpkm == FPKM_THRESHOLD:
    #             # I use a dict to filter becasuse it is faster than a list
    #             dict_uniq_transcrip[transcript_unique_id] = transcript
    #             # uniq_list_transcripts.append(transcript)
    #
    # print '#Num of Transcripts over the threshold ', FPKM_THRESHOLD, ' and distinct'
    # print len(dict_uniq_transcrip)
    #
    #
    # # group the transcript to find the alternative events
    #
    # # I used to group by gene_id, but this become problematic when I take data from diferents readings (use diff id for
    # # the same transcripts when no reference gene was found)
    # # then use TSS to group (transcript start splicing)
    #
    # dict_group_transcripts = {}
    #
    # # for the Novel isoforms, the gene id is not consistent, then I will group isoforms for chromosome and TSS
    # # (transcript start site)
    # for transcript in dict_uniq_transcrip.itervalues():
    #     # TSS id
    #     TSS_id = '-' + transcript.chromosome + '-' + str(transcript.exons[0].start)
    #     # if not dict_group_transcripts.has_key(transcript.gene_id):
    #     if not dict_group_transcripts.has_key(TSS_id):
    #         dict_group_transcripts[TSS_id] = [transcript]
    #     else:
    #         dict_group_transcripts[TSS_id].append(transcript)
    #
    # # extract each transcript of this gene
    #
    # aslocation = open('ASlocation.out', 'w')
    # complete = open('complete_transcripts.fasta', 'w')
    #
    # for TSS_id, list_transcripts in dict_group_transcripts.iteritems():
    #     if len(list_transcripts) > 1:
    #         fetch_events(list_transcripts, aslocation, complete)
