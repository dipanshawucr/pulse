import json
from features.features_helpers import create_paths_for_cell_line
from features.uniprot_transmem import get_transmembrane_region_features
from features.uniprot_ptm import get_postranscriptional_modification_features
from features.uniprot_elm_read import get_uniprot_elm_features
from helpers.normalize_unicode_data import normalize_unicode_data


def feature_extract_cell_line(cell_line, pulse_path, preprocess_input_path, feature_extract_output_path):
    # TODO: Could probably be better if you put every element in FEATURES_SETTINGS into a dict?
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
    uniprot_elm_db_location = normalize_unicode_data(features_settings["F_ELM2.inp"])
    uniprot_elm_read_output_location = feature_extract_output_path + '/elm_read.out'
    get_uniprot_elm_features(uniprot_exon_indices_location, uniprot_elm_db_location, uniprot_elm_read_output_location)
    print "Finished getting eukaryotic linear motif scores."

    ########################
    # DISORDEROME FEATURES #
    ########################

    print "Now getting disorderome features..."
    print "Finished getting disorderome features."
