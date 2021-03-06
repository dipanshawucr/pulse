import os
from transcript import Transcript


def create_paths_for_cell_line(pulse_path, cell_line):
    """
    Creates all necessary paths for preprocess:

    Base folder for cell line:
    - ../output/preprocess/<CELL_LINE_NAME>/

    :param cell_line:
    :param pulse_path:
    :return:
    """

    paths_to_create = [pulse_path + r'/output/preprocess/' + cell_line]
    for path in paths_to_create:
        print path
        if not os.path.exists(path):
            os.makedirs(path)


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


def process_transcripts(transcript_file, dict_of_transcripts, fpkm_threshold):
    """
    For every transcript in list of transcripts, extract those transcripts that are unique and have
    an FPKM equal to fpkm_threshold.

    :param transcript_file:
    :param dict_of_transcripts:
    :param fpkm_threshold:
    :return:
    """
    dictionary_of_unique_transcripts = {}
    list_transcripts = dict_of_transcripts[transcript_file]
    for transcript in list_transcripts:
        exon_ids = ''
        for exon in transcript.exons:
            exon_ids += str(exon.start) + '-' + str(exon.end) + '.'

        transcript_unique_id = transcript.chromosome + '-' + exon_ids

        if transcript_unique_id not in dictionary_of_unique_transcripts and transcript.fpkm == fpkm_threshold:
            dictionary_of_unique_transcripts[transcript_unique_id] = transcript

    print 'Number of transcripts over the threshold ', fpkm_threshold, ' and are distinct:'
    print len(dictionary_of_unique_transcripts)
    return {"list_transcripts": list_transcripts,
            "dictionary_of_unique_transcripts": dictionary_of_unique_transcripts}


def fetch_events(transcripts, as_location, complete, events_location):
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

                        if c1 is not None and c2 is not None and c1[1].id == c2[0].id:
                            event_id = c1[0].id + c1[1].id + c2[0].id + c2[1].id + exon_pairs_i[0].id + exon_pairs_i[
                                1].id
                            if event_id not in events:
                                events.append(event_id)
                                exclusion_intron = c1[1].id

                                # Print to events.fa
                                print >> events_location, '>I.' + transcript_j.id + exclusion_intron + '-' + str(c1[0].length) + '=' + str(
                                    c1[1].length) + '=' + str(c2[1].length)
                                print >> events_location, c1[0].nseq + c1[1].nseq + c2[1].nseq

                                print >> events_location, '>E.' + transcript_i.id + exclusion_intron + '-' + str(c1[0].length) + '=' + str(
                                    c1[1].length) + '=' + str(c2[1].length)
                                print >> events_location, c1[0].nseq + c2[1].nseq

                                # Print to AS_location
                                print >> as_location, transcript_j.chromosome.split('r')[-1] + '\t' + str(
                                    c1[0].start) + '\t' + str(c1[0].end) + '\t' + transcript_j.id + exclusion_intron + \
                                '\t' + 'C1' + '\t' + transcript_j.strand
                                print >> as_location, transcript_j.chromosome.split('r')[-1] + '\t' + str(
                                    c1[1].start) + '\t' + str(c1[1].end) + '\t' + transcript_j.id + exclusion_intron + \
                                '\t' + 'A' + '\t' + transcript_j.strand
                                print >> as_location, transcript_j.chromosome.split('r')[-1] + '\t' + str(
                                    c2[1].start) + '\t' + str(c2[1].end) + '\t' + transcript_j.id + exclusion_intron + \
                                '\t' + 'C2' + '\t' + transcript_j.strand
                                print >> as_location, transcript_i.chromosome.split('r')[-1] + '\t' + str(
                                    c1[0].start) + '\t' + str(c1[0].end) + '\t' + transcript_i.id + exclusion_intron + \
                                '\t' + 'C1' + '\t' + transcript_i.strand
                                print >> as_location, transcript_i.chromosome.split('r')[-1] + '\t' + str(
                                    c1[1].start) + '\t' + str(c1[1].end) + '\t' + transcript_i.id + exclusion_intron + \
                                '\t' + 'A' + '\t' + transcript_i.strand
                                print >> as_location, transcript_i.chromosome.split('r')[-1] + '\t' + str(
                                    c2[1].start) + '\t' + str(c2[1].end) + '\t' + transcript_i.id + exclusion_intron + \
                                '\t' + 'C2' + '\t' + transcript_i.strand

                                # Print to complete
                                print >> complete, '>' + transcript_j.id
                                print >> complete, transcript_j.n_seq()
                                print >> complete, '>' + transcript_i.id
                                print >> complete, transcript_i.n_seq()

