import sys
import requests
import json
import time
import unicodedata


class EnsemblRestClient:

    def __init__(self, server='http://grch37.rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, hdrs=None, params=None):

        # check if we need to rate limit
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            response = requests.get(self.server + endpoint, headers={ "Content-Type" : "application/json"}, timeout=3)
            if not response.ok:
                response.raise_for_status()
            return json.loads(response.content)
        except requests.exceptions.ReadTimeout as e1:
            pass
        except requests.HTTPError as e:
            pass

    def get_gene_of_transcript(self, transcript):
        gene = None
        transcripts = self.perform_rest_action(
            '/xrefs/id/{0}'.format(transcript)
        )
        if transcripts:
            display_id = transcripts[1]['display_id']
            retrieved_gene = self.perform_rest_action(
                '/xrefs/symbol/human/{0}'.format(display_id)
            )
            gene = retrieved_gene[0]['id'].encode('utf-8')
        return gene


def run(transcript):
    client = EnsemblRestClient()
    client.get_gene_of_transcript(transcript)

if __name__ == '__main__':
    symbols_list = ['ENST00000450305',
                    'ENST00000469289',
                    'ENST00000456328',
                    'ENST00000408384',
                    'ENST00000515242',
                    'ENST00000518655',
                    'ENST00000473358']

    if len(sys.argv) == 3:
        symbol = sys.argv[1:]
    else:
        symbol = 'ENST00000456328'

    genes_list = []
    for symbol in symbols_list:
        genes_list.append(run(symbol))