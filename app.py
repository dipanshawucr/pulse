# Main command-line interface for PULSE
import os
import sys
import json
import pickle

from for_preprocess import for_preprocess, cufflinks, samtools

from preprocess import preprocess_helpers, blast
from preprocess.reference_genome import load_and_pickle_reference_genome, load_pickled_reference_genome


try:
    PULSE_PATH = os.environ['PULSE_PATH']
except KeyError:
    print 'Environment variable PULSE_PATH not setup'
    print 'Using working directory'
    PULSE_PATH = os.getcwd()

PREPROCESS_SETTINGS = json.load(open(PULSE_PATH + '/preprocess/preprocess_settings.json', 'r'))

if __name__ == "__main__":
    all_cell_lines = os.listdir(PULSE_PATH + "/input/cell_lines")

    ##########################################################################
    # SAMTOOLS AND CUFFLINKS

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

    ##########################################################################

    proceed_to_preprocessing = raw_input("Press Y to proceed to preprocessing: \n")

    if proceed_to_preprocessing == 'Y':

        # Uncomment below to load and pickle reference genome.
        # print "LOADING AND PICKLING REFERENCE GENOME."
        # load_and_pickle_reference_genome(PULSE_PATH, PREPROCESS_SETTINGS)
        # print "SUCCESSFULLY PICKLED REFERENCE GENOME!"

        # LOAD PICKLED REFERENCE GENOME
        print "LOADING PICKLED REFERENCE GENOME (happens only once for all cell lines)..."
        ref_genome = load_pickled_reference_genome(PULSE_PATH)
        print "PICKLED REFERENCE GENOME LOADED!"

        ##########################################################################

        all_cell_lines_for_preprocess = os.listdir(PULSE_PATH + '/output/for_preprocess')

        # Preprocessing: extract AS events, filter map, then generate indices
        for cell_line in all_cell_lines_for_preprocess:
            preprocess_helpers.create_paths_for_transcript(PULSE_PATH, cell_line)
            print "Preprocessing paths created for: " + cell_line
            print "Beginning alternative splicing extraction step for: " + cell_line
            print "Loading assembled transcripts... \n"
            transcript_file_location = PULSE_PATH + '/output/for_preprocess/' + cell_line + \
                                       '/cufflinks_output/transcripts.gtf'

            dict_of_transcripts = dict()

            dict_of_transcripts[cell_line] = \
                preprocess_helpers.load_assembled_transcripts(transcript_file_location,
                                                              ref_genome)

            print "Assembled transcripts for " + cell_line + " loaded!\n"

            ##########################################################################

            # CREATE A DICTIONARY OF ALL UNIQUE TRANSCRIPTS

            FPKM_THRESHOLD = PREPROCESS_SETTINGS["FPKM_THRESHOLD"]

            processed_transcripts = preprocess_helpers.process_transcripts(
                cell_line,
                dict_of_transcripts,
                FPKM_THRESHOLD
            )

            dictionary_of_unique_transcripts = processed_transcripts["dictionary_of_unique_transcripts"]
            list_transcripts = processed_transcripts["list_transcripts"]

            ##########################################################################

            # USE TSS to group and fetch all AS events

            dict_group_transcripts = {}

            for transcript in dictionary_of_unique_transcripts.itervalues():
                # TSS id
                TSS_id = '-' + transcript.chromosome + '-' + str(transcript.exons[0].start)
                # if not dict_group_transcripts.has_key(transcript.gene_id):
                if TSS_id not in dict_group_transcripts:
                    dict_group_transcripts[TSS_id] = [transcript]
                else:
                    dict_group_transcripts[TSS_id].append(transcript)

            as_location_output = open(PULSE_PATH + '/output/preprocess/' +
                                      cell_line + '/ASlocation.out', 'w')
            complete_output = open(PULSE_PATH + '/output/preprocess/' +
                                   cell_line + '/complete_transcripts.fasta', 'w')
            as_events_output = open(PULSE_PATH + '/output/preprocess/' +
                                    cell_line + '/events.fa', 'w')

            print "Fetching all exclusion/inclusion events... \n"
            for TSS_id, list_transcripts in dict_group_transcripts.iteritems():
                if len(list_transcripts) > 1:
                    print >> as_events_output, \
                        preprocess_helpers.fetch_events(list_transcripts,
                                                        as_location_output,
                                                        complete_output,
                                                        as_events_output)

            as_location_output.close()
            complete_output.close()
            as_events_output.close()

            ##########################################################################

            # BLAST
            uniprot_file_location = PREPROCESS_SETTINGS["UNIPROT_FILE_LOCATION"]
            events_file_location = PULSE_PATH + '/output/preprocess/' + cell_line + '/events.fa'
            output_location = PULSE_PATH + '/output/preprocess/' + cell_line + '/blastx_from_AS_events.out'
            blast.blast_events(uniprot_file_location, events_file_location, output_location)

            ##########################################################################

            # filter map

            ##########################################################################

            # generate index

            ##########################################################################
    else:
        print "Invalid input. Quitting..."
        quit()



            # Feature extraction

            # Machine learning