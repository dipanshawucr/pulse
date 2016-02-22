# The most time-expensive tasks are for_preprocessing
# (cufflinks and samtools) and preprocessing (blast)
# For every cell line, creates a new job for beagle using
# submit_job.py

import os
import sys
from for_preprocess import for_preprocess, cufflinks, samtools
from preprocess.reference_genome import load_and_pickle_reference_genome, load_pickled_reference_genome
from preprocess.main import preprocess_cell_line

try:
    PULSE_PATH = os.environ['PULSE_PATH']
except KeyError:
    print('Environment variable PULSE_PATH not setup')
    print('Using working directory')
    PULSE_PATH = os.getcwd()

if __name__ == "__main__":
    all_cell_lines = os.listdir(PULSE_PATH + "/input/cell_lines")
    print("Found cell lines: ")
    print(all_cell_lines)
    print()
    PREPROCESS_SETTINGS = ''

    for cell_line in all_cell_lines:
        for_preprocess.create_paths_for_cell_line(cell_line)
        print("for_preprocessing paths created for: " + cell_line)
        print("Running samtools for: " + cell_line)
        samtools.run_samtools(PULSE_PATH, cell_line)
        print("Finished samtools for: " + cell_line)
        print("Running cufflinks for: " + cell_line)
        cufflinks.run_cufflinks(PULSE_PATH, cell_line)
        print("Finished cufflinks for: " + cell_line)

    # TODO: Should move pickling and loading of reference genome to preprocess/main.py
    # Uncomment below to load and pickle reference genome.
    print("LOADING AND PICKLING REFERENCE GENOME.")
    load_and_pickle_reference_genome(PULSE_PATH, PREPROCESS_SETTINGS)
    print("SUCCESSFULLY PICKLED REFERENCE GENOME!")
    # print("LOADING PICKLED REFERENCE GENOME (happens only once for all cell lines)...")
    ref_genome = load_pickled_reference_genome(PULSE_PATH)
    # print("PICKLED REFERENCE GENOME LOADED!")

    for cell_line in all_cell_lines_for_preprocess:
            PREPROCESS_OUTPUT_PATH = PULSE_PATH + '/output/preprocess/' + cell_line
            preprocess_cell_line(cell_line, ref_genome, PULSE_PATH, PREPROCESS_OUTPUT_PATH)
    else:
        print "Invalid input. Skipping preprocessing step..."
        pass

