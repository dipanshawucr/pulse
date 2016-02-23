# The most time-expensive tasks are for_preprocessing
# (cufflinks and samtools) and preprocessing (blast)
# For every cell line, this script submits a cufflinks task.

# For dependencies, check for_preprocess/cufflinks

import os
from submit_job import submit_new_job
from preprocess.reference_genome import load_and_pickle_reference_genome, load_pickled_reference_genome


try:
    PULSE_PATH = os.environ['PULSE_PATH']
except KeyError:
    print('Environment variable PULSE_PATH not setup')
    print('Using working directory')
    PULSE_PATH = os.getcwd()

if __name__ == "__main__":
    all_cell_lines = os.listdir(PULSE_PATH + "/input/cell_lines")
    print(all_cell_lines)

    PREPROCESS_SETTINGS = PULSE_PATH + '/preprocess/preprocess_settings.json'

    # Uncomment below to load and pickle reference genome.
    print("LOADING AND PICKLING REFERENCE GENOME.")
    load_and_pickle_reference_genome(PULSE_PATH, PREPROCESS_SETTINGS)
    print("SUCCESSFULLY PICKLED REFERENCE GENOME!")

    print("LOADING PICKLED REFERENCE GENOME (happens only once for all cell lines)...")
    ref_genome = load_pickled_reference_genome(PULSE_PATH)
    print("PICKLED REFERENCE GENOME LOADED!")

    for cell_line in all_cell_lines:
        PREPROCESS_OUTPUT_PATH = PULSE_PATH + '/output/preprocess/' + cell_line
        submit_new_job('preprocess/main.py', ['-c', cell_line, '-r', ref_genome, '-p', PULSE_PATH,
                                              '-o', PREPROCESS_OUTPUT_PATH])
