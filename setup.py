# sets up the root folder for PULSE pipeline

import os


def setup_for_samtools_and_cufflinks():
    paths_to_create = [r'./input/cell_lines', r'./output/for_preprocess']
    for path in paths_to_create:
        if not os.path.exists(path):
            os.makedirs(path)


# PREPROCESS SETUP
def setup_for_preprocess():
    paths_to_create = [r'./output', r'./output/preprocess']
    for path in paths_to_create:
        if not os.path.exists(path):
            os.makedirs(path)

if __name__ == "__main__":
    setup_for_samtools_and_cufflinks()
    print("Make sure you add your raw .bam files for your cell lines"
          " into input/cell_lines!")
    setup_for_preprocess()