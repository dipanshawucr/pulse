import os
import sys
from for_preprocess import for_preprocess, cufflinks, samtools

try:
    PULSE_PATH = os.environ['PULSE_PATH']
except KeyError:
    print 'Environment variable PULSE_PATH not setup'
    print 'Using working directory'
    PULSE_PATH = os.getcwd()

if __name__ == "__main__":
    