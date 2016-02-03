from Bio.Seq import Seq
from generate_index_and_libs_functions import translate_backwards, translate, get_n_seq


def generate1(blastx_file, complete_transcripts_location, filtermap_uniprot_exon_indices_map,
              cdna_output_location, pseq_output_location):
    """

    :param blastx_file:
    :param complete_transcripts_location:
    :param filtermap_uniprot_exon_indices_map:
    :param cdna_output_location:
    :param pseq_output_location:
    :return:
    """
    blastx_output = open(blastx_file, 'r')

    filter_index = {}
    with open(filtermap_uniprot_exon_indices_map, 'r') as input_file:
        for line in input_file:
            data = line.split('\t')
            filter_index[data[0]] = 1

    cdna_output = open(cdna_output_location, 'w')
    pseq_output = open(pseq_output_location, 'w')

    for blastx_hit in blastx_output:

        try:
            (query_id, subject_id, perc_identity, aln_length, mismatch_count,
             gap_open_count, query_start, query_end, subject_start,
             subject_end, e_val, bit_score) = blastx_hit.split('\t')
            if query_id in filter_index:

                n_seq = Seq(get_n_seq(query_id, complete_transcripts_location))

                if int(query_start) > int(query_end):

                    pseq, cdna_seq, s, e = translate_backwards(n_seq, int(query_start), int(query_end))
                    # Save AA seq
                    print >> pseq_output, '>' + str(query_id)
                    print >> pseq_output, pseq
                    # Save CDNA seq
                    print >> cdna_output, '>' + str(query_id)
                    print >> cdna_output, cdna_seq

                else:
                    tn_seq = str(n_seq)
                    pseq, cdna_seq = translate(tn_seq, int(query_start))

                    print >> pseq_output, '>' + str(query_id)
                    print >> pseq_output, pseq
                    # Save CDNA seq
                    print >> cdna_output, '>' + str(query_id)
                    print >> cdna_output, cdna_seq

            else:
                pass
        except Exception, e:
            print e

    pseq_output.close()
    cdna_output.close()
