# Main command-line interface for PULSE
import os

# TODO: Could refactor for_preprocess better
from for_preprocess import for_preprocess, cufflinks, samtools
from preprocess.reference_genome import load_and_pickle_reference_genome, load_pickled_reference_genome
from preprocess.main import preprocess_cell_line
from features.main import feature_extract_cell_line

try:
    PULSE_PATH = os.environ['PULSE_PATH']
except KeyError:
    print 'Environment variable PULSE_PATH not setup'
    print 'Using working directory'
    PULSE_PATH = os.getcwd()

if __name__ == "__main__":
    all_cell_lines = os.listdir(PULSE_PATH + "/input/cell_lines")

    ##########################
    # SAMTOOLS AND CUFFLINKS #
    ##########################

    proceed_to_samtools_cufflinks = raw_input("Press Y to proceed to cufflinks and samtools: ")
    if proceed_to_samtools_cufflinks == 'Y':
        for cell_line in all_cell_lines:
            for_preprocess.create_paths_for_cell_line(cell_line)
            print "for_preprocessing paths created for: " + cell_line
            print "Running samtools for: " + cell_line
            samtools.run_samtools(PULSE_PATH, cell_line)
            print "Finished samtools for: " + cell_line
            print "Running cufflinks for: " + cell_line
            cufflinks.run_cufflinks(PULSE_PATH, cell_line)
            print "Finished cufflinks for: " + cell_line
    else:
        "Invalid input. Skipping samtools and cufflinks step..."
        pass

    proceed_to_preprocessing = raw_input("Press Y to proceed to preprocessing: ")

    #################
    # PREPROCESSING #
    #################

    if proceed_to_preprocessing == 'Y':
        # TODO: Should move pickling and loading of reference genome to preprocess/main.py
        # Uncomment below to load and pickle reference genome.
        # print "LOADING AND PICKLING REFERENCE GENOME."
        # load_and_pickle_reference_genome(PULSE_PATH, PREPROCESS_SETTINGS)
        # print "SUCCESSFULLY PICKLED REFERENCE GENOME!"
        print "LOADING PICKLED REFERENCE GENOME (happens only once for all cell lines)..."
        ref_genome = load_pickled_reference_genome(PULSE_PATH)
        print "PICKLED REFERENCE GENOME LOADED!"
        all_cell_lines_for_preprocess = os.listdir(PULSE_PATH + '/output/for_preprocess')
        for cell_line in all_cell_lines_for_preprocess:
            PREPROCESS_OUTPUT_PATH = PULSE_PATH + '/output/preprocess/' + cell_line
            preprocess_cell_line(cell_line, ref_genome, PULSE_PATH, PREPROCESS_OUTPUT_PATH)
    else:
        print "Invalid input. Skipping preprocessing step..."
        pass

    ######################
    # FEATURE EXTRACTION #
    ######################

    print "Now entering feature extraction."
    proceed_to_feature_extraction = raw_input("Press Y to proceed to feature extraction: ")
    if proceed_to_feature_extraction == "Y":
        all_cell_lines_for_features = os.listdir(PULSE_PATH + '/output/preprocess')
        for cell_line in all_cell_lines_for_features:
            PREPROCESS_OUTPUT_PATH = PULSE_PATH + '/output/preprocess/' + cell_line
            FEATURE_EXTRACT_OUTPUT_PATH = PULSE_PATH + '/output/features/' + cell_line
            feature_extract_cell_line(cell_line, PULSE_PATH, PREPROCESS_OUTPUT_PATH, FEATURE_EXTRACT_OUTPUT_PATH)
    else:
        print "Invalid input. Skipping feature extraction..."
        pass

    ###########
    # ML STEP #
    ###########
    print "Now entering machine learning step."
    proceed_to_machine_learning = raw_input("Press Y to proceed to machine learning step: ")
    if proceed_to_machine_learning == "Y":
        pass
    else:
        print "Invalid input. Skipping machine learning step..."
        pass

    # TODO: Port random forest classifier over to Python
