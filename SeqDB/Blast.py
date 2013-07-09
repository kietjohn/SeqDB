from decimal import *
from Bio.Blast import NCBIXML
from parameters import COLONIZER, NON_COLONIZER

getcontext().prec = 2


class BlastRecord(object):

    def __init__(self, ID, raw_blast_result, blast_object=None):
        self.db_index = ID
        self.blast_result = blast_object

        self.pursue = 0
        self.hits = []

        if not self.blast_result:
            with open(raw_blast_result, 'r') as record:
                self.blast_result = NCBIXML.read(record)

        self.clone = self.blast_result.query

        for align in self.blast_result.alignments:
            for hit in align.hsps:
                genome = align.accession
                organism = align.hit_def
                identity = float(Decimal(hit.identities) / Decimal(hit.align_length))
                self.hits.append((self.db_index, genome, organism, hit.expect,
                                hit.query_start, hit.query_end,
                                hit.sbjct_start, hit.sbjct_end, identity))

        self.colonizer_matches = [item for item in self.hits
                if (item[1] in COLONIZER or item[2].split(' ') in COLONIZER)]
        self.non_colonizer_matches = [item for item in self.hits
                if (item[1] in NON_COLONIZER or item[2].split(' ') in NON_COLONIZER)]
        if len(self.non_colonizer_matches) == 0:
            self.pursue = 1

    def match(self):
        first_match = list(self.hits[0])
        first_match.append(self.pursue)
        return tuple(first_match)

    def __str__(self):
        return ('clone ' + self.clone + ', ' +
                str(len(self.hits)) + ' hits' + ', ' +
                str(len(self.non_colonizer_matches)) + ' non_colonizer')

    def __repr__(self):
        return self.hits


def multiple_blast_records(raw_blast_result):
    with open(raw_blast_result, 'r') as handle:
        blast_records = NCBIXML.parse(handle)
    for item in blast_records:
        pass
