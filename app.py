# Main command-line interface for PULSE
import os
from for_preprocess import cufflinks, samtools

try:
    PULSE_PATH = os.environ['PULSE_PATH']
except KeyError:
    print 'Environment variable PULSE_PATH not setup'
    print 'Using working directory'
    PULSE_PATH = os.getcwd()

if __name__ == "__main__":
    all_cell_lines = os.listdir("./input/cell_lines")

     # First run samtools and cufflinks for all cell lines
    for cell_line in all_cell_lines:
        samtools.create_paths_for_cell_line(cell_line)
        print "for_preprocessing paths created for: " + cell_line
        print "Running samtools for: " + cell_line
        samtools.run_samtools(PULSE_PATH, cell_line)
        print "Finished samtools for: " + cell_line
        print "Running cufflinks for: " + cell_line
        cufflinks.run_cufflinks(PULSE_PATH, cell_line)
        print "Finished cufflinks for: " + cell_line

    for cell_line in all_cell_lines:

        # Preprocessing: extract AS events, filter map, then generate indices


        # Feature extraction

        # Machine learning