# sets up the root folder for PULSE pipeline

import os

# TODO: Download the uniprot using API and insert into pulse/input
# After the uniprot_sprot file is downloaded into input/, need to run:
# formatdb -i ./input/uniprot_sprot.fasta -p T
# creates a hashed version of uniprot_sprot database

# TODO: Download the reference genome using an API and insert into pulse/input/
# TODO: Check for dependencies like samtools and cufflinks and blast, and install if necessary


def setup_for_samtools_and_cufflinks():
    paths_to_create = [r'./input/cell_lines', r'./output', r'./output/for_preprocess']
    for path in paths_to_create:
        if not os.path.exists(path):
            print("Path created: ", path)
            os.makedirs(path)


# PREPROCESS SETUP
def setup_for_preprocess():
    paths_to_create = [r'./output', r'./output/preprocess']
    for path in paths_to_create:
        if not os.path.exists(path):
            print("Path created: ", path)
            os.makedirs(path)


# PREPROCESS SETUP
def setup_for_feature_extraction():
    paths_to_create = [r'./output', r'./output/features']
    for path in paths_to_create:
        if not os.path.exists(path):
            print("Path created: ", path)
            os.makedirs(path)


# ML SETUP
def setup_for_machine_learning():
    paths_to_create = [r'./output', r'./output/machine']
    for path in paths_to_create:
        if not os.path.exists(path):
            print("Path created: ", path)
            os.makedirs(path)


if __name__ == "__main__":
    setup_for_samtools_and_cufflinks()
    print("Make sure you add your raw .bam files for your cell lines"
          " into input/cell_lines!")
    setup_for_preprocess()
    setup_for_feature_extraction()
    setup_for_machine_learning()
