"""
Extract exclusion and inclusion events from transcriptome.
"""

from transcript import Transcript

def load_assembled_trans(filename, cromosomes):
    '''Load transcript from aradobsis file, create a object for each and add them to a list_transcripts
    return:
        List of transcrips objects '''

    datafile = open(filename)
    data = datafile.readlines()
    list_transcripts = []

    for line in data:
        column = line.split('\t')
        exon_id = column[8].split(' ')

        if column[2] == "transcript":
            try:
                # print '>'+itranscript.id+itranscript.cromosome
                # print itranscript.nSeq()
                list_transcripts.append(itranscript)


            except:
                pass

            transcript_info = column[8].split(' ')
            transcript_id = transcript_info[3].strip()
            gene_id = transcript_info[1].strip()
            chromosome = column[0]

            itranscript = Transcript(transcript_id, gene_id, chromosome)


        elif column[2] == "exon":
            exon_id = column[8].split(' ')
            exon_start = int(column[3])
            exon_end = int(column[4])

            if itranscript.id == exon_id[3]:
                itranscript.add_exon([exon_start, exon_end],
                                     cromosomes[itranscript.chromosome][exon_start - 1:exon_end])

            else:
                print exon_id, itranscript.id

    return list_transcripts


def load_sequence_(cromosome):
    '''load the whole nucleotic sequence genome of aradob.
    return a dictionary for each cromosome where the sequence is  a string '''

    filename = '../../Data/Genome_hg19/' + cromosome + '.fa'
    datafile = open(filename)
    data = datafile.readlines()
    chr_sequence = ''
    for line in data:
        if line[0] == '>':
            pass
        else:
            chr_sequence += line.strip()

    datafile.close()
    # print chr_sequence[begin:finish]
    # print len(chr_sequence[begin:finish])

    return chr_sequence.upper()


def fetch_events(transcripts):
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

                        if c1 != None and c2 != None and c1[1].id == c2[0].id:
                            # save_event()
                            event_id = c1[0].id + c1[1].id + c2[0].id + c2[1].id + exon_pairs_i[0].id + exon_pairs_i[
                                1].id
                            if not event_id in events:
                                # print 'lol'
                                # print c1[0].id,c1[1].id,c2[0].id,c2[1].id
                                # print exon_pairs_i[0].id,exon_pairs_i[1].id
                                events.append(event_id)
                                # print transcript_j.id, transcript_j.transcript_info
                                # print transcript_i.id, transcript_i.transcript_info
                                exclusion_intron = c1[1].id
                                print '>I.' + transcript_j.id + exclusion_intron + '-' + str(c1[0].length) + '=' + str(
                                    c1[1].length) + '=' + str(c2[1].length)
                                print c1[0].nseq + c1[1].nseq + c2[1].nseq
                                # print event.info, event.exons_ids
                                # print 'exclusion'
                                print '>E.' + transcript_i.id + exclusion_intron + '-' + str(c1[0].length) + '=' + str(
                                    c1[1].length) + '=' + str(c2[1].length)
                                print c1[0].nseq + c2[1].nseq


def is_conserv(putative_splicing_events):
    for event in putative_splicing_events:
        for e in putative_splicing_events:
            if not event.is_equal(e):
                # print event.info, event.exons_ids
                # print '----'
                # print  e.info, e.exons_ids
                # dump_data(putative_splicing_events)
                return False
    return True


def dump_data(putative_splicing_events):
    print 'DUMPING EVENTS'
    for event in putative_splicing_events:
        print event.info, event.exons_ids


# load the aradobisis genome
cromosomes = {'chr1': '', 'chr2': '', 'chr3': '', 'chr4': '', 'chr5': '', \
              'chr6': '', 'chr7': '', 'chr8': '', 'chr9': '', 'chr10': '', \
              'chr11': '', 'chr12': '', 'chr13': '', 'chr14': '', 'chr15': '', \
              'chr16': '', 'chr17': '', 'chr18': '', 'chr19': '', 'chr20': '', 'chr21': '', 'chr22': '', \
              'chrY': '', 'chrX': '', 'chrM': ''}
for cromo in cromosomes.iterkeys():
    cromosomes[cromo] = load_sequence_(cromo)

# read transcript and convert them to splicing events


list_transcripts = load_assembled_trans('../../Data/Normal/NBS_denovo_transcripts.1.gtf', cromosomes)

print len(list_transcripts)


# generate a list with al genes involved

list_gene_id = []

for transcript in list_transcripts:
    if not transcript.gene_id in list_gene_id:
        list_gene_id.append(transcript.gene_id)
    else:
        pass

# extract each transcript of this gen

for gen in list_gene_id:
    gene_transcripts = []

    for x in list_transcripts:
        if x.gene_id == gen:
            gene_transcripts.append(x)

    if len(gene_transcripts) > 1:
        fetch_events(gene_transcripts)
