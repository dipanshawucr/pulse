import json
import os
import time
from preprocess import preprocess_helpers, blast, filtermap, generate_index_and_libs


# TODO: Can also break this into smaller bits
def preprocess_cell_line(cell_line, ref_genome, pulse_path, preprocess_output_path):
    # TODO: This alternative splicing event detection can be made simpler
    PREPROCESS_SETTINGS = json.load(open(pulse_path + '/preprocess/preprocess_settings.json', 'r'))
    preprocess_helpers.create_paths_for_cell_line(pulse_path, cell_line)
    print "Preprocessing paths created for: " + cell_line
    print "Beginning alternative splicing extraction step for: " + cell_line
    print "Loading assembled transcripts... \n"
    transcript_file_location = pulse_path + '/output/for_preprocess/' + cell_line + \
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

    as_location_output = open(preprocess_output_path + '/as_location.out', 'w')
    complete_output = open(preprocess_output_path + '/complete_transcripts.fasta', 'w')
    as_events_output = open(preprocess_output_path + '/events.fa', 'w')
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

    # ##########################################################################
    # BLAST
    # TODO: Automate the step where formatdb has to be called
    uniprot_file_location = PREPROCESS_SETTINGS["UNIPROT_FILE_LOCATION"]
    events_file_location = pulse_path + '/output/preprocess/' + cell_line + '/events.fa'
    output_location = pulse_path + '/output/preprocess/' + cell_line + '/blastx_from_as_events.out'
    print "Now blasting the events.fa file from: " + cell_line + "\n"
    blast_exit_code = blast.blast_events(uniprot_file_location, events_file_location, output_location)
    print blast_exit_code
    ##########################################################################

    blast_exit_code = 0
    if blast_exit_code == 0:
        print "Blast success!"

        # Generate mapping between uniprot to splicing events with filter map
        blast_file = output_location
        uniprot_fasta = uniprot_file_location
        isoform_fasta = preprocess_output_path + '/complete_transcripts.fasta'

        if not os.path.exists(preprocess_output_path + '/temp'):
            os.makedirs(preprocess_output_path + '/temp')

        filtermap_not_len_collapsed_output = preprocess_output_path + \
                                             '/temp/rnaseq_huniprot_corrected_len.txt'
        filtermap_len_collapsed_output = preprocess_output_path + \
                                         '/temp/rnaseq_huniprot_corrected_len_collapsed.txt'
        filtermap_uniprot_exon_indices_map = preprocess_output_path + '/uniprot_exon_indices_map.out'

        print "Running filtermap 1\n"
        filtermap.filter_map1(blast_file, filtermap_not_len_collapsed_output)
        time.sleep(1)
        print "Running filtermap 2\n"
        filtermap.filter_map2(filtermap_not_len_collapsed_output, filtermap_len_collapsed_output)
        time.sleep(1)
        print "Running filtermap 3\n"
        uniprot_ddbb = filtermap.filter_map3(uniprot_fasta)
        time.sleep(1)
        print "Running filtermap 4\n"
        isoforms_ddbb = filtermap.filter_map4(isoform_fasta)
        time.sleep(1)
        print "Running filtermap 5\n"
        has_mapping = filtermap.filter_map5(filtermap_len_collapsed_output)
        time.sleep(3)
        print "Running filtermap 6\n"
        filtermap.filter_map6(filtermap_not_len_collapsed_output, filtermap_uniprot_exon_indices_map,
                              has_mapping, uniprot_ddbb, isoforms_ddbb)

        ##########################################################################
        # generate index
        blastx_file = preprocess_output_path + '/temp/rnaseq_huniprot_corrected_len_collapsed.txt'
        complete_transcripts = isoform_fasta
        cdna_output_location = preprocess_output_path + '/cdna_transcripts.fasta'
        pseq_output_location = preprocess_output_path + '/p_seq_isoforms.fas'
        print "Running generate 1"
        generate_index_and_libs.generate1(blastx_file, complete_transcripts,
                                          filtermap_uniprot_exon_indices_map,
                                          cdna_output_location, pseq_output_location)

        ##########################################################################
    else:
        print "BLAST Failed."
        # TODO: Write better error statement for debugging.
