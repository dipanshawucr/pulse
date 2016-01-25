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
