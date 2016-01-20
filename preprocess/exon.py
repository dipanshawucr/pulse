class Exon(object):

    def __init__(self, transcript, chr_map, sequence=''):
        self.nseq = sequence
        self.start = int(chr_map[0])
        self.end = int(chr_map[1])
        self.id = transcript.chromosome + '_' + str(self.start) + '_' + str(self.end)
        self.length = self.end - self.start + 1
