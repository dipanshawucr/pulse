import json
import time
from features.features_helpers import create_paths_for_cell_line
from features.uniprot_transmem import get_transmembrane_region_features
from features.uniprot_ptm import get_postranscriptional_modification_features
from features.uniprot_elm_read import get_uniprot_elm_features
from features.generate_iupred_file import generate_iupred_file
from features.uniprot_disorder import get_uniprot_disorder_features
from features.uniprot_domain_read import get_uniprot_domain_read
from features.run_pfam_scan import start_pfam_scan
from features.uniprot_core import get_sable_scores
from features.mutation_features import get_mutation_features
from features.conservation_conversion_query import create_query_file
from features.convert_from_hg19_hg18 import use_remap_api
from features.seq_conserv import generate_sequence_conservation_features
from features.event_conserv import generate_event_conservation_feature_table
from features.generate_ml_input import generate_machine_learning_matrix
from features.network_features import generate_network_features
from helpers.normalize_unicode_data import normalize_unicode_data


def feature_extract_cell_line(cell_line, pulse_path, preprocess_input_path, feature_extract_output_path):

    # TODO: Could probably be better if you put every element in FEATURES_SETTINGS into a dict?
    # TODO: Should move all of param_files_locations to be relative to pulse!
    features_settings = json.load(open(pulse_path + '/features/features_settings.json'))
    create_paths_for_cell_line(pulse_path, cell_line)
    print "Features paths created for: " + cell_line

    #########################
    # TRANSMEMBRANE SCORING #
    #########################

    print "Now getting transmembrane region features..."
    uniprot_exon_indices_location = preprocess_input_path + '/uniprot_exon_indices_map.out'

    uniprot_tm_indices_db_location = normalize_unicode_data(features_settings["F_UNIPROT_TRANSMEM_INDICES_LOCATION"])
    uniprot_tm_read_output_location = feature_extract_output_path + '/transmem_read.out'
    get_transmembrane_region_features(uniprot_exon_indices_location, uniprot_tm_indices_db_location,
                                      uniprot_tm_read_output_location)
    print "Finished getting transmembrane region features."
    time.sleep(2)
    ################
    # PTM FEATURES #
    ################

    print "Now getting post-transcriptional modifications..."
    uniprot_ptm_db_location = normalize_unicode_data(features_settings["F_PTMS_LOCATION"])
    uniprot_ptm_read_output_location = feature_extract_output_path + '/ptm_read.out'
    get_postranscriptional_modification_features(uniprot_exon_indices_location, uniprot_ptm_db_location,
                                                 uniprot_ptm_read_output_location)
    print "Finished getting post-transcriptional modification features."

    ###################################
    # EUKARYOTIC LINEAR MOTIF SCORING #
    ###################################

    print "Now getting eukaryotic linear motif scores..."
    uniprot_elm_db_location = normalize_unicode_data(features_settings["F_ELM2_LOCATION"])
    uniprot_elm_read_output_location = feature_extract_output_path + '/elm_read.out'
    get_uniprot_elm_features(uniprot_exon_indices_location, uniprot_elm_db_location, uniprot_elm_read_output_location)
    print "Finished getting eukaryotic linear motif scores."

    #############################
    # DISORDEROME HELPER SCRIPT #
    #############################

    print "Now running helper files for disorderome..."
    p_seq_output_location = preprocess_input_path + '/p_seq_isoforms.fas'
    iupred_isoforms_output_location = feature_extract_output_path + '/iupred_isoforms.out'
    iupred_install_path = normalize_unicode_data(features_settings["IUPRED_INSTALL_PATH"])
    generate_iupred_file(p_seq_output_location, feature_extract_output_path,
                          iupred_install_path, iupred_isoforms_output_location)
    print "Now done running helper file for disorderome."

    ########################
    # DISORDEROME FEATURES #
    ########################

    print "Now getting disorderome features..."
    canonical_db_location = pulse_path + '/input/info_canonical_v3.ddbb'
    disorder_read_out_location = feature_extract_output_path + '/disorder_read.out'
    get_uniprot_disorder_features(pulse_path, uniprot_exon_indices_location, iupred_isoforms_output_location,
                                  canonical_db_location, disorder_read_out_location)
    print "Finished getting disorderome features."

    ##########################
    # PFAM & DOMAIN FEATURES #
    ##########################

    print "Now running pfam_scan..."
    pfam_scan_script_location = pulse_path + '/helpers/pfam_scan.pl'
    pfam_input_location = preprocess_input_path + '/p_seq_isoforms.fas'
    pfam_output_location = feature_extract_output_path + '/pfam_done.out'
    hmmer3_data_location = normalize_unicode_data(features_settings["HMMER3_DATA_LOCATION"])
    pfam_exit_code = start_pfam_scan(pfam_scan_script_location, pfam_input_location,
                                     pfam_output_location, hmmer3_data_location)
    pfam_exit_code = 0
    if pfam_exit_code == 0:
        print "pfam_scan successful"
        print "Now getting uniprot domain features..."
        f_pfam_special_db_location = normalize_unicode_data(features_settings["F_PFAM_SPECIAL_LOCATION"])
        domain_read_output_location = feature_extract_output_path + '/domain_read.out'
        get_uniprot_domain_read(f_pfam_special_db_location, canonical_db_location, uniprot_exon_indices_location,
                                pfam_output_location, domain_read_output_location)
        print "Finished getting uniprot domain features."

        #################
        # SABLE SCORING #
        #################

        print "Now getting features for SABLE..."
        f_sable_db_location = normalize_unicode_data(features_settings["F_SABLE_LOCATION"])
        uniprot_core_output_location = feature_extract_output_path + '/core_read.out'
        get_sable_scores(uniprot_exon_indices_location, f_sable_db_location, uniprot_core_output_location)
        print "Finished getting features for SABLE."

        ####################
        # MUTATION SCORING #
        ####################

        print "Now getting mutation scores..."
        f_mutations_db_location = normalize_unicode_data(features_settings["F_MUTATIONS_LOCATION"])
        mutation_features_output_location = feature_extract_output_path + '/mutation_read.out'
        get_mutation_features(uniprot_exon_indices_location, f_mutations_db_location,
                              mutation_features_output_location)
        print "Finished getting mutation scores."

        #################################
        # CONSERVATION/NETWORK FEATURES #
        #################################
        print "Creating query file conservation/network features..."
        conservation_query_output_location = feature_extract_output_path + '/conservation_query.txt'
        as_location_file = preprocess_input_path + '/as_location.out'
        create_query_file(as_location_file, conservation_query_output_location)
        print "Finished creating query file."

        print "Converting between hg19 to hg18 using the query file generated..."
        remap_api_location = pulse_path + '/helpers/remap_api.pl'
        remap_mode = "asm-asm"
        remap_from = "GCF_000001405.13"
        remap_to = "GCF_000001405.12"
        remap_input_location = conservation_query_output_location
        remap_output_location = feature_extract_output_path + '/report_conservationQuery.txt.xls'
        remap_exit_code = use_remap_api(remap_api_location, remap_mode, remap_from, remap_to,
                                        remap_input_location, remap_output_location)
        print "Conversion complete."

        remap_exit_code = 0
        if remap_exit_code == 0:
            print "Now generating sequence conservation features..."
            f_phastcons_db_location = normalize_unicode_data(features_settings["F_PHASTCONS_HG18_BED_LOCATION"])
            as_location_file = preprocess_input_path + '/as_location.out'
            remapped_coordinates_file = feature_extract_output_path + '/report_conservationQuery.txt'
            sequence_conservation_output_location = feature_extract_output_path + '/sequenceCon_read.out'
            generate_sequence_conservation_features(f_phastcons_db_location, as_location_file,
                                                    remapped_coordinates_file, sequence_conservation_output_location)
            print "Finished generating sequence conservation features."

            print "Now generating event_conservation feature table..."
            f_event_conservation_db_location = normalize_unicode_data(features_settings["F_EVENT_CONSERVATION_LOCATION"])
            event_conservation_output = feature_extract_output_path + '/eventCon_read.out'
            generate_event_conservation_feature_table(f_phastcons_db_location, f_event_conservation_db_location,
                                                      as_location_file, remapped_coordinates_file,
                                                      event_conservation_output)
            print "Finished generating event_conservation feature table."

            print "Now generating network features using gene names..."
            f_uniprot_genewiki_location = normalize_unicode_data(features_settings["F_UNIPROT_GENEWIKI_LOCATION"])
            f_degree_location = normalize_unicode_data(features_settings["F_DEGREE_LOCATION"])
            network_features_output_location = feature_extract_output_path + '/degree_read.out'
            generate_network_features(f_uniprot_genewiki_location, f_degree_location, uniprot_exon_indices_location,
                                      network_features_output_location)
            print "Finished generating network features."

            ###############################
            # GENERATE MATRIX OF FEATURES #
            ###############################

            print "Now generating ML input..."
            machine_learning_output_path = pulse_path + '/output/machine/' + cell_line
            generated_ml_output = machine_learning_output_path + '/features_headers.txt'
            generated_ml_names = machine_learning_output_path + '/names.txt'
            ts_read = normalize_unicode_data(features_settings["ML_AS_EVENT_CLASSIFICATION"])
            generate_machine_learning_matrix(uniprot_exon_indices_location, uniprot_core_output_location,
                                             network_features_output_location,
                                             disorder_read_out_location, domain_read_output_location,
                                             uniprot_elm_read_output_location, event_conservation_output,
                                             mutation_features_output_location, uniprot_ptm_read_output_location,
                                             sequence_conservation_output_location, uniprot_tm_read_output_location,
                                             ts_read, generated_ml_output, generated_ml_names)
            print "Finished generating ML input."

        else:
            print "Remapping failed"
            exit()
    else:
        print "pfam_scan failed"
        exit()
