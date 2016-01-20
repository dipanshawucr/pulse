from exon import Exon

class Transcript(object):

    def __init__(self, trans_id, gene_id, chro):

        # self.transcript_info = parsedline[8].split(' ')
        self.id = trans_id.strip()
        self.gene_id = gene_id.strip()
        self.chromosome = chro.strip()
        self.exons = []

    def add_exon(self, chr_map, sequence=''):

        exon = Exon(self, chr_map, sequence)

        self.exons.append(exon)

        return

    def nSeq(self):

        nseq = ''
        for e in self.exons:
            nseq += e.nseq

        return nseq

    def extract_events_splicing(self):

        events = []

        n = len(self.exons)

        i = 0
        while True:
            list_exons = self.exons[i:i + 3]
            i += 1
            if len(list_exons) == 3:
                event = Splicing_Event(list_exons, self.transcript_info)
                events.append(event)
            else:
                break

        return events

    def extract_exon_ids(self):

        self.list_exon_ids = []

        for e in self.exons:
            self.list_exon_ids.append(e.id)

        return

    def iterpairs(self):

        for x in range(len(self.exons) - 1):
            pair = [self.exons[x], self.exons[x + 1]]
            if len(pair) == 2:
                yield pair
            else:
                return