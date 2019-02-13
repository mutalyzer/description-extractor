import extractor
import pytest


REFERENCE = 'ACGTCGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT'


TESTS = [
    ('(empty)',
     '', '',
     []),

    ('=',
     REFERENCE, REFERENCE,
     [{'source': 'reference', 'location': {'end': {'position': 44, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}]),

    ('7A>G',
     REFERENCE, 'ACGTCGGTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{'source': 'reference', 'location': {'end': {'position': 6, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 7, 'type': 'point'}, 'type': 'range', 'start': {'position': 6, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 7, 'type': 'point'}, 'type': 'range', 'start': {'position': 6, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 44, 'type': 'point'}, 'type': 'range', 'start': {'position': 7, 'type': 'point'}}, 'type': 'equal'}]),

    ('7del',
     REFERENCE, 'ACGTCGTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{'source': 'reference', 'location': {'end': {'position': 6, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 6, 'type': 'point'}, 'type': 'range', 'start': {'position': 6, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 7, 'type': 'point'}, 'type': 'range', 'start': {'position': 6, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 44, 'type': 'point'}, 'type': 'range', 'start': {'position': 7, 'type': 'point'}}, 'type': 'equal'}]),

    ('7_8del',
     REFERENCE, 'ACGTCGTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{'source': 'reference', 'location': {'end': {'position': 6, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 6, 'type': 'point'}, 'type': 'range', 'start': {'position': 6, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 8, 'type': 'point'}, 'type': 'range', 'start': {'position': 6, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 44, 'type': 'point'}, 'type': 'range', 'start': {'position': 8, 'type': 'point'}}, 'type': 'equal'}]),

    ('6_7insC',
     REFERENCE, 'ACGTCGCATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{'source': 'reference', 'location': {'end': {'position': 6, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 7, 'type': 'point'}, 'type': 'range', 'start': {'position': 6, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 6, 'type': 'point'}, 'type': 'range', 'start': {'position': 6, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 44, 'type': 'point'}, 'type': 'range', 'start': {'position': 6, 'type': 'point'}}, 'type': 'equal'}]),

    ('6_7insCC',
     REFERENCE, 'ACGTCGCCATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{'source': 'reference', 'location': {'end': {'position': 6, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 8, 'type': 'point'}, 'type': 'range', 'start': {'position': 6, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 6, 'type': 'point'}, 'type': 'range', 'start': {'position': 6, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 44, 'type': 'point'}, 'type': 'range', 'start': {'position': 6, 'type': 'point'}}, 'type': 'equal'}]),

    ('7_11inv',
     REFERENCE, 'ACGTCGCGAATCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{'source': 'reference', 'location': {'end': {'position': 6, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'source': 'reference', 'location': {'end': {'position': 11, 'type': 'point'}, 'type': 'range', 'start': {'position': 6, 'type': 'point'}}, 'type': 'inversion'}, {'source': 'reference', 'location': {'end': {'position': 44, 'type': 'point'}, 'type': 'range', 'start': {'position': 11, 'type': 'point'}}, 'type': 'equal'}]),

    ('6_7inv',
     REFERENCE, 'ACGTCTCTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{'source': 'reference', 'location': {'end': {'position': 5, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'source': 'reference', 'location': {'end': {'position': 7, 'type': 'point'}, 'type': 'range', 'start': {'position': 5, 'type': 'point'}}, 'type': 'inversion'}, {'source': 'reference', 'location': {'end': {'position': 44, 'type': 'point'}, 'type': 'range', 'start': {'position': 7, 'type': 'point'}}, 'type': 'equal'}]),

    ('7delinsCC',
     REFERENCE, 'ACGTCGCCTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{'source': 'reference', 'location': {'end': {'position': 6, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 8, 'type': 'point'}, 'type': 'range', 'start': {'position': 6, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 7, 'type': 'point'}, 'type': 'range', 'start': {'position': 6, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 44, 'type': 'point'}, 'type': 'range', 'start': {'position': 7, 'type': 'point'}}, 'type': 'equal'}]),

    ('21_23delinsTTTT',
     REFERENCE, 'ACGTCGATTCGCTAGCTTCGTTTTGATAGATAGAGATATAGAGAT',
     [{'source': 'reference', 'location': {'end': {'position': 20, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 24, 'type': 'point'}, 'type': 'range', 'start': {'position': 20, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 23, 'type': 'point'}, 'type': 'range', 'start': {'position': 20, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 44, 'type': 'point'}, 'type': 'range', 'start': {'position': 23, 'type': 'point'}}, 'type': 'equal'}]),

    ('7_8delinsTC',
     REFERENCE, 'ACGTCGTCTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{'source': 'reference', 'location': {'end': {'position': 6, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 8, 'type': 'point'}, 'type': 'range', 'start': {'position': 6, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 8, 'type': 'point'}, 'type': 'range', 'start': {'position': 6, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 44, 'type': 'point'}, 'type': 'range', 'start': {'position': 8, 'type': 'point'}}, 'type': 'equal'}]),

    ('7dup',
     REFERENCE, 'ACGTCGAATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{'source': 'reference', 'location': {'end': {'position': 7, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 8, 'type': 'point'}, 'type': 'range', 'start': {'position': 7, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 7, 'type': 'point'}, 'type': 'range', 'start': {'position': 7, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 44, 'type': 'point'}, 'type': 'range', 'start': {'position': 7, 'type': 'point'}}, 'type': 'equal'}]),

    ('6_7dup',
     REFERENCE, 'ACGTCGAGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{'source': 'reference', 'location': {'end': {'position': 7, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 9, 'type': 'point'}, 'type': 'range', 'start': {'position': 7, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 7, 'type': 'point'}, 'type': 'range', 'start': {'position': 7, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 44, 'type': 'point'}, 'type': 'range', 'start': {'position': 7, 'type': 'point'}}, 'type': 'equal'}]),

    ('5_7dup',
     REFERENCE, 'ACGTCGACGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
     [{'source': 'reference', 'location': {'end': {'position': 7, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 10, 'type': 'point'}, 'type': 'range', 'start': {'position': 7, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 7, 'type': 'point'}, 'type': 'range', 'start': {'position': 7, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 44, 'type': 'point'}, 'type': 'range', 'start': {'position': 7, 'type': 'point'}}, 'type': 'equal'}]),

    ('[5_6insTT;17del;26A>C;35dup]',
     'ATGATGATCAGATACAGTGTGATACAGGTAGTTAGACAA', 'ATGATTTGATCAGATACATGTGATACCGGTAGTTAGGACAA',
     [{'source': 'reference', 'location': {'end': {'position': 5, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 7, 'type': 'point'}, 'type': 'range', 'start': {'position': 5, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 5, 'type': 'point'}, 'type': 'range', 'start': {'position': 5, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 16, 'type': 'point'}, 'type': 'range', 'start': {'position': 5, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 18, 'type': 'point'}, 'type': 'range', 'start': {'position': 18, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 17, 'type': 'point'}, 'type': 'range', 'start': {'position': 16, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 25, 'type': 'point'}, 'type': 'range', 'start': {'position': 17, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 27, 'type': 'point'}, 'type': 'range', 'start': {'position': 26, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 26, 'type': 'point'}, 'type': 'range', 'start': {'position': 25, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 34, 'type': 'point'}, 'type': 'range', 'start': {'position': 26, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 36, 'type': 'point'}, 'type': 'range', 'start': {'position': 35, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 34, 'type': 'point'}, 'type': 'range', 'start': {'position': 34, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 39, 'type': 'point'}, 'type': 'range', 'start': {'position': 34, 'type': 'point'}}, 'type': 'equal'}]),

    ('[26A>C;30C>A;35G>C]',
     'TAAGCACCAGGAGTCCATGAAGAAGATGGCTCCTGCCATGGAATCCCCTACTCTACTGTG', 'TAAGCACCAGGAGTCCATGAAGAAGCTGGATCCTCCCATGGAATCCCCTACTCTACTGTG',
     [{'source': 'reference', 'location': {'end': {'position': 25, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 26, 'type': 'point'}, 'type': 'range', 'start': {'position': 25, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 26, 'type': 'point'}, 'type': 'range', 'start': {'position': 25, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 29, 'type': 'point'}, 'type': 'range', 'start': {'position': 26, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 30, 'type': 'point'}, 'type': 'range', 'start': {'position': 29, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 30, 'type': 'point'}, 'type': 'range', 'start': {'position': 29, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 34, 'type': 'point'}, 'type': 'range', 'start': {'position': 30, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 35, 'type': 'point'}, 'type': 'range', 'start': {'position': 34, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 35, 'type': 'point'}, 'type': 'range', 'start': {'position': 34, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 60, 'type': 'point'}, 'type': 'range', 'start': {'position': 35, 'type': 'point'}}, 'type': 'equal'}]),

    ('[26_29inv;30C>G]',
     'TAAGCACCAGGAGTCCATGAAGAAGATGGCTCCTGCCATGGAATCCCCTACTCTA', 'TAAGCACCAGGAGTCCATGAAGAAGCCATGTCCTGCCATGGAATCCCCTACTCTA',
     [{'source': 'reference', 'location': {'end': {'position': 25, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'source': 'reference', 'location': {'end': {'position': 29, 'type': 'point'}, 'type': 'range', 'start': {'position': 25, 'type': 'point'}}, 'type': 'inversion'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 30, 'type': 'point'}, 'type': 'range', 'start': {'position': 29, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 30, 'type': 'point'}, 'type': 'range', 'start': {'position': 29, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 55, 'type': 'point'}, 'type': 'range', 'start': {'position': 30, 'type': 'point'}}, 'type': 'equal'}]),

    ('[26_29inv;30C>G;41del]',
     'TAAGCACCAGGAGTCCATGAAGAAGATGGCTCCTGCCATGGAATCCCCTACTCTA', 'TAAGCACCAGGAGTCCATGAAGAAGCCATGTCCTGCCATGAATCCCCTACTCTA',
     [{'source': 'reference', 'location': {'end': {'position': 25, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'source': 'reference', 'location': {'end': {'position': 29, 'type': 'point'}, 'type': 'range', 'start': {'position': 25, 'type': 'point'}}, 'type': 'inversion'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 30, 'type': 'point'}, 'type': 'range', 'start': {'position': 29, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 30, 'type': 'point'}, 'type': 'range', 'start': {'position': 29, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 39, 'type': 'point'}, 'type': 'range', 'start': {'position': 30, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 39, 'type': 'point'}, 'type': 'range', 'start': {'position': 39, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 40, 'type': 'point'}, 'type': 'range', 'start': {'position': 39, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 55, 'type': 'point'}, 'type': 'range', 'start': {'position': 40, 'type': 'point'}}, 'type': 'equal'}]),

    ('37_38ins15_24',
     'ATGGCGGCGGTGGTCGCCCTCTCCTTGAGGCGCCGGTTGCCGGCCACAACCCTTGGCGGA', 'ATGGCGGCGGTGGTCGCCCTCTCCTTGAGGCGCCGGTCGCCCTCTCCTGCCGGCCACAACCCTTGGCGGA',
     [{'source': 'reference', 'location': {'end': {'position': 37, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'reference', 'location': {'end': {'position': 24, 'type': 'point'}, 'type': 'range', 'start': {'position': 14, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 37, 'type': 'point'}, 'type': 'range', 'start': {'position': 37, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 60, 'type': 'point'}, 'type': 'range', 'start': {'position': 37, 'type': 'point'}}, 'type': 'equal'}]),

    ('37_38ins3_26inv',
     'ATGGCGGCGGTGGTCGCCCTCTCCTTGAGGCGCCGGTTGCCGGCCACAACCCTTGGCGGA', 'ATGGCGGCGGTGGTCGCCCTCTCCTTGAGGCGCCGGTAAGGAGAGGGCGACCACCGCCGCCTGCCGGCCACAACCCTTGGCGGA',
     [{'source': 'reference', 'location': {'end': {'position': 37, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'reference', 'inverted': True, 'location': {'end': {'position': 26, 'type': 'point'}, 'type': 'range', 'start': {'position': 2, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 37, 'type': 'point'}, 'type': 'range', 'start': {'position': 37, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 60, 'type': 'point'}, 'type': 'range', 'start': {'position': 37, 'type': 'point'}}, 'type': 'equal'}]),

    ('18_42dup',
     'ATGGCGGCGGTGGTCGCCCTCTCCTTGAGGCGCCGGTTGCCGGCCACAACCCTTGGCGGA', 'ATGGCGGCGGTGGTCGCCCTCTCCTTGAGGCGCCGGTTGCCGCCTCTCCTTGAGGCGCCGGTTGCCGGCCACAACCCTTGGCGGA',
     [{'source': 'reference', 'location': {'end': {'position': 42, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'reference', 'location': {'end': {'position': 42, 'type': 'point'}, 'type': 'range', 'start': {'position': 17, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 42, 'type': 'point'}, 'type': 'range', 'start': {'position': 42, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 60, 'type': 'point'}, 'type': 'range', 'start': {'position': 42, 'type': 'point'}}, 'type': 'equal'}]),

    ('42_43ins[CA;19_42]',
     'ATGGCGGCGGTGGTCGCCCTCTCCTTGAGGCGCCGGTTGCCGGCCACAACCCTTGGCGGA', 'ATGGCGGCGGTGGTCGCCCTCTCCTTGAGGCGCCGGTTGCCGCACTCTCCTTGAGGCGCCGGTTGCCGGCCACAACCCTTGGCGGA',
     [{'source': 'reference', 'location': {'end': {'position': 42, 'type': 'point'}, 'type': 'range', 'start': {'position': 0, 'type': 'point'}}, 'type': 'equal'}, {'inserted': [{'source': 'observed', 'location': {'end': {'position': 44, 'type': 'point'}, 'type': 'range', 'start': {'position': 42, 'type': 'point'}}}, {'source': 'reference', 'location': {'end': {'position': 42, 'type': 'point'}, 'type': 'range', 'start': {'position': 18, 'type': 'point'}}}], 'source': 'reference', 'location': {'end': {'position': 42, 'type': 'point'}, 'type': 'range', 'start': {'position': 42, 'type': 'point'}}, 'type': 'deletion_insertion'}, {'source': 'reference', 'location': {'end': {'position': 60, 'type': 'point'}, 'type': 'range', 'start': {'position': 42, 'type': 'point'}}, 'type': 'equal'}]),
]


@pytest.mark.parametrize('_, reference, observed, variants', TESTS)
def test_variants(_, reference, observed, variants):
    assert variants == extractor.describe_dna(reference, observed)
