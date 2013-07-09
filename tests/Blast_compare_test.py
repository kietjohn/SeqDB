from nose.tools import *
from SeqDB.Blast import BlastRecord
from Bio.Blast import NCBIXML

def blast_test():
    '''BLAST result interpretation

    Given several BLAST result xml, load the results, and then do the sorting.
    Compare the generated result with manual input result.
    '''
    blast_object1 = BlastRecord(55, 'tests/test_data/blast/single_blast1.xml')
    assert_equal(blast_object1.match(), (55, 'HM991502', 'Pseudomonas fluorescens strain Q8r1-96 type III secretion gene cluster, complete sequence', 0.0, 20, 992, 11854, 10866, 0.98, 1))

    blast_object2 = BlastRecord(46, 'tests/test_data/blast/single_blast2.xml')
    assert_equal(blast_object2.match(), (46, 'CP002585', 'Pseudomonas brassicacearum subsp. brassicacearum NFM421, complete genome', 0.0, 19, 525, 636636, 636111, 0.96, 0))

    with open('tests/test_data/blast/single_blast2.xml', 'r') as handle:
        blast = NCBIXML.read(handle)
    multi_test = BlastRecord(65, 'dummy_place_holder', blast)
    assert_equal(multi_test.match(), (65, 'CP002585', 'Pseudomonas brassicacearum subsp. brassicacearum NFM421, complete genome', 0.0, 19, 525, 636636, 636111, 0.96, 0))


