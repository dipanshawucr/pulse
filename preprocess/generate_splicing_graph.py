import pickle
from splicing_graph import SplicingGraph, SplicingSite, Edge
import pprint
from ensembl_client import EnsemblRestClient

cuff1 = open('cuff1test.txt', 'r')


def find_lowest_exon_position(transcripts):
    """
    Returns the lowest exon position out of all given exons for a transcript.

    :param exons:
    :return:
    """
    all_positions = set()
    for transcript in transcripts:
        for exon in transcripts[transcript]['exons']:
            for exon_pos in transcripts[transcript]['exons'][exon]:
                all_positions.add(exon_pos)
    return min(all_positions)


def clean_cufflinks_file(cufflinks_file):
    """
    This method generates the full dictionary of genes to transcripts to exons (including strand sense)
    from the raw cufflinks output.

    :param cufflinks_file:
    :return:
    """
    full_dictionary = {}
    for line in cufflinks_file:
        tokens = line.split('\t')
        gene_id = tokens[8].split()[1].translate(None, '";')
        transcript_id = tokens[8].split()[3].translate(None, '";')
        start = int(tokens[3])
        end = int(tokens[4])
        sense = tokens[6]
        if tokens[2] == 'transcript':
            if gene_id not in full_dictionary:
                full_dictionary[gene_id] = {transcript_id: {'sense': sense, 'exons': {}}}
            else:
                full_dictionary[gene_id][transcript_id] = {'sense': sense, 'exons': {}}
        elif tokens[2] == 'exon':
            exon_number = int(tokens[8].split()[5].translate(None, '";'))
            full_dictionary[gene_id][transcript_id]['exons'][exon_number] = (start, end)
    return full_dictionary


def search_ensembl_from_transcript_to_gene(raw_cufflinks_genes_to_transcript_to_exons, ensembl_client):
    """
    For each transcript, check if they have an entry in Ensembl.
    Check - does transcript have a gene match in Ensembl?
    Yes -> Add to gene_to_transcripts_to_exons under its retrieved gene.
    No -> Adds to gene_to_transcripts_to_exons under None.

    :param raw_cufflinks_genes_to_transcript_to_exons:
    :return:
    """

    genes_to_transcripts_to_exons = {}
    for gene in raw_cufflinks_genes_to_transcript_to_exons:
        for transcript in raw_cufflinks_genes_to_transcript_to_exons[gene]:
            print("Now processing transcript: ", transcript)
            ensembl_gene = ensembl_client.get_gene_of_transcript(transcript)
            if ensembl_gene not in genes_to_transcripts_to_exons:
                genes_to_transcripts_to_exons[ensembl_gene] = \
                    [{transcript: raw_cufflinks_genes_to_transcript_to_exons[gene][transcript]}]
            else:
                genes_to_transcripts_to_exons[ensembl_gene] += \
                    [{transcript: raw_cufflinks_genes_to_transcript_to_exons[gene][transcript]}]
    return genes_to_transcripts_to_exons


def get_all_splice_sites(genes_to_transcripts_to_exons_dict):
    """
    Tracks all of the splice sites for each gene.

    :param genes_to_transcripts_to_exons_dict:
    :return:
    """
    gene_to_splice_site = {}
    for gene in genes_to_transcripts_to_exons_dict:
        print(gene)
        if gene not in gene_to_splice_site:
            gene_to_splice_site[gene] = []
        else:
            pass

        for transcript in genes_to_transcripts_to_exons[gene]:
            transcript = dict(transcript)
            for exon in transcript.values():
                for tuple in exon['exons'].values():
                    if tuple[0] not in gene_to_splice_site[gene]:
                        gene_to_splice_site[gene] += [tuple[0]]
                    if tuple[1] not in gene_to_splice_site[gene]:
                        gene_to_splice_site[gene] += [tuple[1]]

        gene_to_splice_site[gene].sort()
        print("Number of splice sites for " + str(gene) + ": " + str(len(gene_to_splice_site[gene])))

    return gene_to_splice_site

def build_splicing_graph_for_gene(sg, gene, gene_to_transcripts_to_exons):
    """
    Takes all of a gene's transcripts as input and creates a splicing graph object with nodes as
    all of the splice sites and edges labelled as either introns or exons.

    :return:
    """
    # add root node
    sg.add_node(SplicingSite(None, None, 'root'))

    # build graph
    for transcript in genes_to_transcripts_to_exons[gene]:
        transcript_id = transcript.keys()[0]
        # print(transcript_id)
        # print(transcript_id, transcript[transcript_id]['exons'])
        for exon in transcript[transcript_id]['exons']:
            splice_sites = transcript[transcript_id]['exons'][exon]
            node1 = SplicingSite(splice_sites[0], transcript_id, 'start')
            sg.add_node(node1)
            node2 = SplicingSite(splice_sites[1], transcript_id, 'end')
            sg.add_node(node2)
            sg.add_edge(Edge(True, node1, node2))

    # add leaf node
    sg.add_node(SplicingSite(None, None, 'leaf'))

    # create edges for:
    # the introns,
    # between last node and artificial leaf node,
    # and between first node and artificial root

    return sg


if __name__ == '__main__':

    raw_cufflinks_dictionary = clean_cufflinks_file(cuff1)

    # Uncomment this block to re-pickle!
    # client = EnsemblRestClient()
    # genes_to_transcripts_to_exons = search_ensembl_from_transcript_to_gene(raw_cufflinks_dictionary, client)
    #
    # # pickle the result of genes_to_transcripts_to_exons to save time for testing
    # with open('test_mapping_of_genes_to_transcripts_to_exons.pkl', 'wb') as output:
    #     pickle.dump(genes_to_transcripts_to_exons, output, pickle.HIGHEST_PROTOCOL)

    # Load pickled data
    with open('test_mapping_of_genes_to_transcripts_to_exons.pkl', 'rb') as file_to_load:
        genes_to_transcripts_to_exons = pickle.load(file_to_load)

    # # TODO could probably just get the start site of the gene...
    # for gene in raw_cufflinks_dictionary:
    #     last_seen_splice_site = find_lowest_exon_position(raw_cufflinks_dictionary[gene])
    #     print(last_seen_splice_site)

    pp = pprint.PrettyPrinter(indent=1)
    # pp.pprint(genes_to_transcripts_to_exons)

    # get all splice sites
    print(get_all_splice_sites(genes_to_transcripts_to_exons))

    # Create objects
    gene_graphs = {key: SplicingGraph(key) for (key, value) in genes_to_transcripts_to_exons.iteritems()}

    for gene in gene_graphs:
        # print("Current gene graph being built: " + str(gene))
        gene_graphs[gene] = build_splicing_graph_for_gene(gene_graphs[gene], gene, genes_to_transcripts_to_exons)
        # print("Gene graph after building: " + str(gene_graphs[gene]) + "\n")

    # print(gene_graphs['ENSG00000223972'])

    # # for all of the nodes in the edges that are starts and ends, check if there exists any overlap
    # dict_of_same_end_sites = {}
    # dict_of_same_start_sites = {}
    # for node in all_splicing_graphs['ENSG00000223972'].nodes():
    #     if node.site_type == 'start':
    #         if node.position not in dict_of_same_start_sites:
    #             dict_of_same_start_sites[node.position] = node.transcript_ids
    #         else:
    #             dict_of_same_start_sites[node.position] += node.transcript_ids
    #     elif node.site_type == 'end':
    #         if node.position not in dict_of_same_start_sites:
    #             dict_of_same_end_sites[node.position] = node.transcript_ids
    #         else:
    #             dict_of_same_end_sites[node.position] += node.transcript_ids
    #
    # print("Same start sites: ")
    # print(dict_of_same_start_sites)
    # print("                     ")
    # print("Same end sites: ")
    # print(dict_of_same_end_sites)
    #
    # print("-----------------")
    # print(all_splicing_graphs['ENSG00000223972'])
    # # Start creating the nodes for the Splicing Graph
    # for gene in genes_to_transcripts_to_exons:
    #     # build_splicing_graph_for_gene(gene)
    #     print(gene, genes_to_transcripts_to_exons[gene])