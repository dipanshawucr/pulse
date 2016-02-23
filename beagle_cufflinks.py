# The most time-expensive tasks are for_preprocessing
# (cufflinks and samtools) and preprocessing (blast)
# For every cell line, this script submits a cufflinks task.

# For dependencies, check for_preprocess/cufflinks

import os
from submit_job import submit_new_job
from for_preprocess import for_preprocess, cufflinks, samtools


try:
    PULSE_PATH = os.environ['PULSE_PATH']
except KeyError:
    print('Environment variable PULSE_PATH not setup')
    print('Using working directory')
    PULSE_PATH = os.getcwd()

if __name__ == "__main__":
    all_cell_lines = os.listdir(PULSE_PATH + "/input/cell_lines")
    print(all_cell_lines)
    for cell_line in all_cell_lines:
        for_preprocess.create_paths_for_cell_line(cell_line)
        submit_new_job(PULSE_PATH + '/for_preprocess/cufflinks.py', ['-p', PULSE_PATH, '-c', cell_line])
