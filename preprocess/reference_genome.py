import pickle
from preprocess.preprocess_helpers import normalize_unicode_data, load_sequence


class ReferenceGenome(object):

    def __init__(self):
        self.chromosomes_dict = {}
        self.chromosomes = {'1': '', '2': '', '3': '', '4': '', '5': '',
                            '6': '', '7': '', '8': '', '9': '', '10': '',
                            '11': '', '12': '', '13': '', '14': '', '15': '',
                            '16': '', '17': '', '18': '', '19': '', '20': '',
                            '21': '', '22': '', 'Y': '', 'X': '', 'MT': ''}

    def add_genome(self, chromosome, sequence):
        self.chromosomes_dict[chromosome] = sequence


def load_and_pickle_reference_genome(pulse_path, preprocess_settings):
    # LOAD REFERENCE GENOME
    ref_genome = ReferenceGenome()
    REFERENCE_GENOME_DIRECTORY = normalize_unicode_data(
        preprocess_settings["REFERENCE_GENOME_DIRECTORY"]
    )

    for chromosome in ref_genome.chromosomes:
        sequence = load_sequence(REFERENCE_GENOME_DIRECTORY, chromosome)
        ref_genome.add_genome(chromosome, sequence)

    # PICKLE REFERENCE GENOME
    with open(pulse_path + '/input/GRCh37_refGenome_pickle.pkl', 'wb') as output:
        pickle.dump(ref_genome, output, pickle.HIGHEST_PROTOCOL)


def load_pickled_reference_genome(pulse_path):
    with open(pulse_path + '/input/GRCh37_refGenome_pickle.pkl', 'rb') as object_to_load:
        return pickle.load(object_to_load)
