# Main command-line interface for PULSE
import os
from auto_samtools_cufflinks import run_samtools, run_cufflinks

try:
    PULSE_PATH = os.environ['PULSE_PATH']
except KeyError:
    print 'Environment variable PULSE_PATH not setup'
    print 'Using working directory'
    PULSE_PATH = os.getcwd()

# First run samtools and cufflinks

# Preprocessing: extract AS events, filter map, then generate indices

# Feature extraction

# Machine learning