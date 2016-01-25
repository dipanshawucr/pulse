import unicodedata
from transcript import Transcript


def normalize_unicode_data(data):
    """
    Returns a Python-normalized String object from unicode.

    :param data:
    :return:
    """
    normalized_data = unicodedata.normalize('NFKD', data).encode('ascii', 'ignore')
    return normalized_data


def load_sequence(refgenome_directory, chromosome):
    filename = refgenome_directory + chromosome + '.fa'
    datafile = open(filename)
    data = datafile.readlines()
    chr_sequence = ''
    for line in data:
        if line[0] == '>':
            pass
        else:
            chr_sequence += line.strip()
    datafile.close()
    return chr_sequence.upper()


def load_assembled_transcripts(filename, ref_genome):
    datafile = open(filename)
    data = datafile.readlines()
    list_transcripts = []
    for line in data:
        column = line.split('\t')
        if column[2] == "transcript":
            try:
                list_transcripts.append(itranscript)
            except:
                pass
            transcript_info = column[8].split(' ')
            transcript_id = transcript_info[3].strip()
            gene_id = transcript_info[1].strip()
            chromosome = column[0]
            if chromosome in ref_genome.chromosomes_dict:
                fpkm = transcript_info[7].strip(";").strip("\"")
                sign = column[6].strip()
                itranscript = Transcript(transcript_id, gene_id, chromosome, fpkm, sign)
            else:
                pass
        elif column[2] == "exon":
            if column[0] in ref_genome.chromosomes_dict:
                transcript_info = column[8].split(' ')
                exon_start = int(column[3])
                exon_end = int(column[4])
                if itranscript.id == transcript_info[3]:
                    itranscript.add_exon([exon_start, exon_end],
                                         ref_genome.chromosomes_dict[itranscript.chromosome][exon_start - 1:exon_end])
                else:
                    print 'WARNING ', transcript_info, itranscript.id
            else:
                pass
    return list_transcripts


def fetch_events(transcripts, aslocation, complete):
    events = []

    for transcript_i in transcripts:
        for exon_pairs_i in transcript_i.iterpairs():
            for transcript_j in transcripts:
                if transcript_i.id != transcript_j.id:
                    c1 = None
                    c2 = None
                    exon_pair = False

                    for exon_pairs_j in transcript_j.iterpairs():
                        if exon_pairs_i[0].id == exon_pairs_j[0].id and exon_pairs_i[1].id == exon_pairs_j[1].id:
                            exon_pair = True
                        else:
                            pass

                    if not exon_pair:
                        for exon_pairs_j in transcript_j.iterpairs():
                            if exon_pairs_i[0].id == exon_pairs_j[0].id:
                                c1 = exon_pairs_j
                            elif exon_pairs_i[1].id == exon_pairs_j[1].id:
                                c2 = exon_pairs_j
                            else:
                                pass

                        if c1 is None and c2 is None and c1[1].id == c2[0].id:
                            event_id = c1[0].id + c1[1].id + c2[0].id + c2[1].id + exon_pairs_i[0].id + exon_pairs_i[
                                1].id
                            if event_id not in events:
                                events.append(event_id)
                                exclusion_intron = c1[1].id
                                print '>I.' + transcript_j.id + exclusion_intron + '-' + str(c1[0].length) + '=' + str(
                                    c1[1].length) + '=' + str(c2[1].length)
                                print c1[0].nseq + c1[1].nseq + c2[1].nseq

                                print '>E.' + transcript_i.id + exclusion_intron + '-' + str(c1[0].length) + '=' + str(
                                    c1[1].length) + '=' + str(c2[1].length)
                                print c1[0].nseq + c2[1].nseq
                                print >> aslocation, transcript_j.chromosome.split('r')[-1] + '\t' + str(
                                    c1[0].start) + '\t' + str(c1[
                                                                  0].end) + '\t' + transcript_j.id + exclusion_intron + '\t' + 'C1' + '\t' + transcript_j.strand
                                print >> aslocation, transcript_j.chromosome.split('r')[-1] + '\t' + str(
                                    c1[1].start) + '\t' + str(c1[
                                                                  1].end) + '\t' + transcript_j.id + exclusion_intron + '\t' + 'A' + '\t' + transcript_j.strand
                                print >> aslocation, transcript_j.chromosome.split('r')[-1] + '\t' + str(
                                    c2[1].start) + '\t' + str(c2[
                                                                  1].end) + '\t' + transcript_j.id + exclusion_intron + '\t' + 'C2' + '\t' + transcript_j.strand
                                print >> aslocation, transcript_i.chromosome.split('r')[-1] + '\t' + str(
                                    c1[0].start) + '\t' + str(c1[
                                                                  0].end) + '\t' + transcript_i.id + exclusion_intron + '\t' + 'C1' + '\t' + transcript_i.strand
                                print >> aslocation, transcript_i.chromosome.split('r')[-1] + '\t' + str(
                                    c1[1].start) + '\t' + str(c1[
                                                                  1].end) + '\t' + transcript_i.id + exclusion_intron + '\t' + 'A' + '\t' + transcript_i.strand
                                print >> aslocation, transcript_i.chromosome.split('r')[-1] + '\t' + str(
                                    c2[1].start) + '\t' + str(c2[
                                                                  1].end) + '\t' + transcript_i.id + exclusion_intron + '\t' + 'C2' + '\t' + transcript_i.strand
                                print >> complete, '>' + transcript_j.id
                                print >> complete, transcript_j.nSeq()
                                print >> complete, '>' + transcript_i.id
                                print >> complete, transcript_i.nSeq()
