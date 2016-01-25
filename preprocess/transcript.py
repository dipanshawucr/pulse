from exon import Exon


class Transcript(object):

    def __init__(self, trans_id, gene_id, chro, fpkm=0.0, sign='+'):

        # self.transcript_info = parsedline[8].split(' ')
        self.id = trans_id.strip()
        self.gene_id = gene_id.strip()
        self.chromosome = chro.strip()
        self.exons = []
        self.fpkm = float(fpkm)
        self.strand = sign

    def add_exon(self, chr_map, sequence=''):
        exon = Exon(self, chr_map, sequence)
        self.exons.append(exon)
        return

    def n_seq(self):
        nucleotide_seq = ''
        for exon in self.exons:
            nucleotide_seq += exon.nseq
        return nucleotide_seq

    def extract_exon_ids(self):
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