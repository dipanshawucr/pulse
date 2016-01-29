# Main command-line interface for PULSE
import os
import sys
import json
import pickle

from for_preprocess import for_preprocess, cufflinks, samtools

from preprocess import preprocess_helpers
from preprocess.reference_genome import ReferenceGenome, load_and_pickle_reference_genome, load_pickled_reference_genome


try:
    PULSE_PATH = os.environ['PULSE_PATH']
except KeyError:
    print 'Environment variable PULSE_PATH not setup'
    print 'Using working directory'
    PULSE_PATH = os.getcwd()

PREPROCESS_SETTINGS = json.load(open('./preprocess/preprocess_settings.json', 'r'))

if __name__ == "__main__":
    # all_cell_lines = os.listdir("./input/cell_lines")
    #
    # # First run samtools and cufflinks for all cell lines
    # for cell_line in all_cell_lines:
    #     for_preprocess.create_paths_for_cell_line(cell_line)
    #     print "for_preprocessing paths created for: " + cell_line
    #     print "Running samtools for: " + cell_line
    #     samtools.run_samtools(PULSE_PATH, cell_line)
    #     print "Finished samtools for: " + cell_line
    #     print "Running cufflinks for: " + cell_line
    #     cufflinks.run_cufflinks(PULSE_PATH, cell_line)
    #     print "Finished cufflinks for: " + cell_line
    #
    # print "Finished cufflinks and samtools for all cell lines: " + str(all_cell_lines)

    proceed_to_preprocessing = raw_input("Press Y to proceed to preprocessing: \n")

    if proceed_to_preprocessing == 'Y':
        ##########################################################################

        # Uncomment to load and pickle reference genome.
        print "LOADING AND PICKLING REFERENCE GENOME."
        load_and_pickle_reference_genome(PULSE_PATH, PREPROCESS_SETTINGS)
        print "SUCCESSFULLY PICKLED REFERENCE GENOME!"

        # LOAD PICKLED REFERENCE GENOME
        print "LOADING PICKLED REFERENCE GENOME (happens once only for all cell lines)..."
        refGenome = load_pickled_reference_genome(PULSE_PATH)
        print "PICKLED REFERENCE GENOME LOADED!"

        ##########################################################################

        all_cell_lines_for_preprocess = os.listdir('./output/for_preprocess')

        # Preprocessing: extract AS events, filter map, then generate indices
        for cell_line in all_cell_lines_for_preprocess:
            preprocess_helpers.create_paths_for_transcript(PULSE_PATH, cell_line)
            print "Preprocessing paths created for: " + cell_line
            print "Beginning alternative splicing extraction step for: " + cell_line
            print "Loading assembled transcripts... \n"
            transcript_file_location = './output/for_preprocess/' + cell_line + '/transcripts.gtf'

            dict_of_transcripts = dict()
            dict_of_transcripts[cell_line] = preprocess_helpers.load_assembled_transcripts(transcript_file_location, ref_genome)
            print "Assembled transcripts for " + cell_line + " loaded!\n"

            print dict_of_transcripts[cell_line]
            # CREATE A DICTIONARY OF ALL UNIQUE TRANSCRIPTS

            # FPKM_THRESHOLD = PREPROCESS_SETTINGS["FPKM_THRESHOLD"]
            #
            # processed_transcripts = preprocess_helpers.process_transcripts(
            #     list_of_transcript_files,
            #     dict_of_transcripts,
            #     FPKM_THRESHOLD
            # )
            #
            # dictionary_of_unique_transcripts = processed_transcripts["dictionary_of_unique_transcripts"]
            # list_transcripts = processed_transcripts["list_transcripts"]

            #
            # # group the transcript to find the alternative events
            #
            # # I used to group by gene_id, but this become problematic
            # # when I take data from different readings (use diff id for
            # # the same transcripts when no reference gene was found)
            # # then use TSS to group (transcript start splicing)
            #
            # dict_group_transcripts = {}
            #
            # # for the Novel isoforms, the gene id is not consistent, then I will group
            # # isoforms for chromosome and TSS (transcript start site)
            #
            # for transcript in dictionary_of_unique_transcripts.itervalues():
            #     # TSS id
            #     TSS_id = '-' + transcript.chromosome + '-' + str(transcript.exons[0].start)
            #     # if not dict_group_transcripts.has_key(transcript.gene_id):
            #     if TSS_id not in dict_group_transcripts:
            #         dict_group_transcripts[TSS_id] = [transcript]
            #     else:
            #         dict_group_transcripts[TSS_id].append(transcript)
            #
            # # extract each transcript of this gene
            #
            # as_location_output = open(normalize_unicode_data(PREPROCESS_SETTINGS["AS_LOCATION_OUTPUT_LOCATION"]) +
            #                   'ASlocation.out', 'w')
            # complete_output = open(normalize_unicode_data(PREPROCESS_SETTINGS["LOCATION_OF_COMPLETE_TRANSCRIPTS"]) +
            #                        'complete_transcripts.fasta', 'w')
            # as_events_output = open(normalize_unicode_data(PREPROCESS_SETTINGS["LOCATION_OF_EVENTS"]) +
            #                         'events.fa', 'w')
            #
            # for TSS_id, list_transcripts in dict_group_transcripts.iteritems():
            #     if len(list_transcripts) > 1:
            #         print >> as_events_output, fetch_events(list_transcripts,
            #                                                 as_location_output,
            #                                                 complete_output)
            #
            # as_location_output.close()
            # complete_output.close()
            # as_events_output.close()
    else:
        print "Invalid input. Quitting..."
        quit()



            # Feature extraction

            # Machine learning